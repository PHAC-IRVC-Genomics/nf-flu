include { fluPrefix } from './functions'

process MINIMAP2 {
  tag "$sample|$segment|$ref_id"
  label 'process_low'

  conda 'bioconda::minimap2=2.28 bioconda::samtools=1.20'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/mulled-v2-058de387f9917a7a63953f496cdd203bca83b790:86215829f86df9201683956877a19d025261ff66-0'
  } else {
    container 'quay.io/biocontainers/mulled-v2-058de387f9917a7a63953f496cdd203bca83b790:86215829f86df9201683956877a19d025261ff66-0'
  }

  input:
  tuple val(sample), val(segment), val(ref_id), path(ref_fasta), path(reads)

  output:
  tuple val(sample), val(segment), val(ref_id), path(ref_fasta), path('*.{bam,bam.bai}'), emit: alignment
  path '*.{flagstat,idxstats,stats}', emit: stats
  path('*.minimap2.log'), emit: log
  path "versions.yml" , emit: versions

  script:
  def prefix        = fluPrefix(sample, segment, ref_id)
  bam               = "${prefix}.bam"
  raw_bam           = "${prefix}.unfiltered.bam"
  flagstat          = "${prefix}.flagstat"
  idxstats          = "${prefix}.idxstats"
  stats             = "${prefix}.stats"
  minimap2_log      = "${prefix}.minimap2.log"

  // Mapping preset based on platform
  def map_option = params.platform == 'nanopore' ? 'map-ont' : 'sr'

  // One combined -F flag (default 2308 = 4 + 256 + 2048: unmapped + secondary + supplementary)
  def flagF = (params.samtools_flagF == null) ? 2308 : params.samtools_flagF as int

  // Threshold for aligned (M, X and =) query bases
  def min_len = (params.min_map_len == null) ? 0 : params.min_map_len as int

  // DVG filtering: Maximum allowed skip/deletion size in CIGAR (N or D operations)
  // DVG reads typically have N operations of 1000s of bp
  // Set to null or 0 to disable this filter
  def max_skip = (params.max_skip_size != null) ? params.max_skip_size as int : 200
  
  // DVG filtering: limit gap size during alignment (default: no limit)
  // Set to 500-1000 to prevent DVG reads with large internal deletions from aligning
  // Set to null or 0 to disable
  def max_gap = (params.max_gap_size != null && params.max_gap_size > 0) ? params.max_gap_size as int : 0
  
  // Disable secondary alignments to ensure single best mapping per read
  boolean allow_secondary = false
  if (params.allow_secondary instanceof Boolean) {
    allow_secondary = params.allow_secondary
  } else if (params.allow_secondary != null) {
    allow_secondary = params.allow_secondary.toString().toLowerCase() in ['1','true','yes','y']
  }

  boolean filter_bam = false
  if (params.filter_bam instanceof Boolean) {
    filter_bam = params.filter_bam
  } else if (params.filter_bam != null) {
    filter_bam = params.filter_bam.toString().toLowerCase() in ['1','true','yes','y']
  }

  // AWK counts only M, X and = CIGAR ops as aligned bases: ignores soft clips, insertions, deletions and skips.
  // Also filters reads with large N/D operations (DVG reads) and ensures single contiguous alignment blocks.
  // Updated to be mawk-compatible with explicit integer conversion and character-by-character parsing.

  """
  MAP_OPTION=$map_option
  FLAG_F=$flagF
  MIN_LEN=$min_len
  MAX_SKIP=$max_skip
  MAX_GAP=$max_gap
  ALLOW_SECONDARY=$allow_secondary
  FILTER_BAM=$filter_bam

  if [ "${params.platform}" = "nanopore" ] && [ "\$MIN_LEN" -gt 0 ] && [ "\$FILTER_BAM" = "true" ]; then
    echo "DVG filtering ENABLED: platform=${params.platform}, MIN_LEN=$min_len, MAX_SKIP=$max_skip, MAX_GAP=$max_gap, ALLOW_SECONDARY=$allow_secondary" >&2
    
    # Build minimap2 options
    MM2_OPTS="-ax \$MAP_OPTION -t${task.cpus}"
    if [ "\$MAX_GAP" -gt 0 ]; then
      MM2_OPTS="\$MM2_OPTS -G \$MAX_GAP"
    fi
    if [ "\$ALLOW_SECONDARY" = "false" ]; then
      MM2_OPTS="\$MM2_OPTS --secondary=no"
    fi
    
    minimap2 \$MM2_OPTS $ref_fasta $reads \\
      | samtools view -h -@${task.cpus} -F \$FLAG_F \\
      | tee >(samtools view -@${task.cpus} -b \\
              | samtools sort -@${task.cpus} -O BAM -o $raw_bam) \\
      | awk -v min_len=\$MIN_LEN -v max_skip=\$MAX_SKIP 'BEGIN {OFS="\t"} {
          if (\$0 ~ /^@/) {print; next}
          
          aligned_length = 0
          max_deletion = 0
          num_blocks = 0
          in_block = 0
          cigar = \$6
          num_str = ""
          
          # Parse CIGAR string character by character
          for (i = 1; i <= length(cigar); i++) {
              c = substr(cigar, i, 1)
              if (c ~ /[0-9]/) {
                  num_str = num_str c
              } else {
                  if (num_str != "") {
                      num_val = int(num_str)
                      
                      # Count M, =, X as aligned bases
                      if (c == "M" || c == "=" || c == "X") {
                          aligned_length += num_val
                          # Track alignment blocks
                          if (in_block == 0) {
                              num_blocks++
                              in_block = 1
                          }
                      }
                      # Track large deletions/skips (DVG signature)
                      else if (c == "N" || c == "D") {
                          if (num_val > max_deletion) {
                              max_deletion = num_val
                          }
                          # Large gap breaks the alignment block
                          if (num_val > 50) {
                              in_block = 0
                          }
                      }
                      
                      num_str = ""
                  }
              }
          }
          
          # Filter criteria:
          # 1. Sufficient aligned length (>= min_len)
          # 2. No large deletions/skips (< max_skip) - filters DVGs (if max_skip > 0)
          # 3. Single contiguous alignment block - filters split alignments
          if (aligned_length >= min_len && (max_skip == 0 || max_deletion <= max_skip) && num_blocks == 1) {
              print
          }
        }' \\
      | samtools view -@${task.cpus} -b \\
      | samtools sort -@${task.cpus} -O BAM > $bam

    samtools index $raw_bam

  else
    echo "AWK-based mapped read-length filtering DISABLED: (platform=${params.platform}, MIN_LEN=\$MIN_LEN)" >&2
    
    # Build minimap2 options even when filtering is disabled
    MM2_OPTS="-ax \$MAP_OPTION -t${task.cpus}"
    if [ "\$MAX_GAP" -gt 0 ]; then
      MM2_OPTS="\$MM2_OPTS -G \$MAX_GAP"
    fi
    if [ "\$ALLOW_SECONDARY" = "false" ]; then
      MM2_OPTS="\$MM2_OPTS --secondary=no"
    fi
    
    minimap2 \$MM2_OPTS $ref_fasta $reads \\
      | samtools view -h -@${task.cpus} -F \$FLAG_F \\
      | samtools view -@${task.cpus} -b \\
      | samtools sort -@${task.cpus} -O BAM > $bam
  fi

  samtools index $bam
  samtools stats $bam > $stats
  samtools flagstat $bam > $flagstat
  samtools idxstats $bam > $idxstats

  ln -s .command.log $minimap2_log

  ${params.platform == 'illumina' ? "samtools faidx $ref_fasta" : ""}

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      minimap2: \$(minimap2 --version 2>&1)
      samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
  END_VERSIONS
  """
}
