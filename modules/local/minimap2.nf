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

  boolean filter_bam = false
  if (params.filter_bam instanceof Boolean) {
    filter_bam = params.filter_bam
  } else if (params.filter_bam != null) {
    filter_bam = params.filter_bam.toString().toLowerCase() in ['1','true','yes','y']
  }

  // AWK counts only M, X and = CIGAR ops as aligned bases: ignores soft clips, insertions, deletions and skips.

  """
  MAP_OPTION=$map_option
  FLAG_F=$flagF
  MIN_LEN=$min_len
  FILTER_BAM=$filter_bam

  if [ "${params.platform}" = "nanopore" ] && [ "\$MIN_LEN" -gt 0 ] && [ "\$FILTER_BAM" = "true" ]; then
    echo "AWK-based mapped read-length filtering ENABLED: (platform=${params.platform}, MIN_LEN=$min_len, filter_bam=$filter_bam)" >&2
    minimap2 -ax \$MAP_OPTION -t${task.cpus} $ref_fasta $reads \\
      | samtools view -h -@${task.cpus} -F \$FLAG_F \\
      | tee >(samtools view -@${task.cpus} -b \\
              | samtools sort -@${task.cpus} -O BAM -o $raw_bam) \\
      | awk 'BEGIN {OFS="\t"} {
          if (\$0 ~ /^@/) {print; next}
          cigar=\$6; aligned_length=0;
          match(cigar, /([0-9]+[M=X])/)
          while (RLENGTH > -1) {
              match_value = substr(cigar, RSTART, RLENGTH - 1) # -1 as we can chop off the cigar operator
              aligned_length += match_value;
              cigar = substr(cigar, RSTART + RLENGTH);
              match(cigar, /([0-9]+[M=])/)
          }
          if (aligned_length >= $min_len) print
        }' \\
      | samtools view -@${task.cpus} -b \\
      | samtools sort -@${task.cpus} -O BAM > $bam

    samtools index $raw_bam

  else
    echo "AWK-based mapped read-length filtering DISABLED: (platform=${params.platform}, MIN_LEN=\$MIN_LEN)" >&2
    minimap2 -ax \$MAP_OPTION -t${task.cpus} $ref_fasta $reads \\
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