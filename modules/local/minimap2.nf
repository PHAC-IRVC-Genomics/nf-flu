include { fluPrefix } from './functions'

process MINIMAP2 {
  tag "$sample|$segment|$ref_id"
  label 'process_low'

  conda 'bioconda::minimap2=2.28 bioconda::samtools=1.20'
  container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? 'https://depot.galaxyproject.org/singularity/mulled-v2-058de387f9917a7a63953f496cdd203bca83b790:86215829f86df9201683956877a19d025261ff66-0' : 'quay.io/biocontainers/mulled-v2-058de387f9917a7a63953f496cdd203bca83b790:86215829f86df9201683956877a19d025261ff66-0' }"

  input:
  tuple val(sample), val(segment), val(ref_id), path(ref_fasta), path(reads)

  output:
  tuple val(sample), val(segment), val(ref_id), path(ref_fasta), path("${bam}"), path("${bam}.bai"), emit: alignment
  tuple val(sample), val(segment), val(ref_id), path(ref_fasta), path("${raw_bam}"), path("${raw_bam}.bai"), emit: raw_alignment, optional: true
  path '*.{flagstat,idxstats,stats}', emit: stats
  path('*.minimap2.log'), emit: log
  path('*.dvg_filter_stats.json'), optional: true, emit: dvg_stats
  path "versions.yml" , emit: versions

  script:
  def prefix        = fluPrefix(sample, segment, ref_id)
  bam               = "${prefix}.bam"
  raw_bam           = "${prefix}.unfiltered.bam"
  flagstat          = "${prefix}.flagstat"
  idxstats          = "${prefix}.idxstats"
  stats             = "${prefix}.stats"
  minimap2_log      = "${prefix}.minimap2.log"
  dvg_stats         = "${prefix}.dvg_filter_stats.json"

  // Mapping preset based on platform
  def map_option = params.platform == 'nanopore' ? 'map-ont' : 'sr'

  // Raw/unfiltered BAM flag: default 4 (unmapped only) — retains secondaries,
  // supplementary, and all other alignments to keep the BAM as raw as possible.
  // Override with params.samtools_flagF_raw if needed.
  // Note: if output_unmapped_reads=true, this flag is set to 0 (no filtering).
  def flagF_raw      = (params.samtools_flagF_raw      != null) ? params.samtools_flagF_raw      as int : 4

  // Filtered BAM flag: 2308 = 4 (unmapped) + 256 (secondary) + 2048 (supplementary)
  // Applied to the filtered BAM only.
  // Override with params.samtools_flagF_filtered if needed.
  def flagF_filtered = (params.samtools_flagF_filtered != null) ? params.samtools_flagF_filtered as int : 2308

  // Normalise segment name once; used by both lookup tables below.
  // Strips leading numeric prefix from IRMA-style names (e.g. '1_PB2' -> 'PB2', '4_HA' -> 'HA').
  // Guards against null values and ensures consistent key matching.
  def segmentUpper = segment?.toString()?.toUpperCase()?.replaceAll(/^\d+_/, '') ?: ''

  // Minimum aligned (M, X and =) query bases required to retain a read in the
  // filtered BAM. Acts as a mapping confidence / read quality floor — independent
  // of DVG filtering. Segment-aware defaults target ~50% segment coverage, which
  // ensures a read contributes meaningfully to consensus/variant calling without
  // being so strict that it discards legitimate reads on short segments.
  //   PB2/PB1     (~2341 bp) -> 1170 bp
  //   PA          (~2151 bp) -> 1075 bp
  //   HA          (~1778 bp) ->  890 bp
  //   NP          (~1565 bp) ->  780 bp
  //   NA          (~1413 bp) ->  700 bp
  //   MP          (~1027 bp) ->  510 bp
  //   NS          ( ~890 bp) ->  445 bp
  // Override globally with params.min_map_len, or accept per-segment defaults.
  // null = use segment-aware defaults; 0 = disabled.
  def minlen_defaults = [
    'PB2': 1170, 'PB1': 1170, 'PA': 1075,
    'HA':   890, 'NP':   780,
    'NA':   700,
    'MP':   510, 'NS':   445
  ]
  def min_len = (params.min_map_len != null)
      ? params.min_map_len as int
      : minlen_defaults.getOrDefault(segmentUpper, 445)

  // Maximum soft-clip ratio allowed in the filtered BAM (0.0–1.0).
  // Reads where soft_clip_bases / total_query_bases exceeds this are discarded.
  // Targets internally-mapping DVG-origin reads that show no large CIGAR gap
  // but have heavily clipped flanks. 0.0 = disabled.
  def max_clip_frac = (params.max_clip_frac != null) ? params.max_clip_frac as double : 0.15

  // DVG CIGAR filter: maximum allowed single N or D operation in the filtered BAM.
  // Thresholds are set at ~20% of segment length sitting above the largest known
  // legitimate biological deletions (~120 bp for NA stalk deletions), while remaining
  // below canonical DVG sizes.
  //   PB2/PB1     (~2341 bp) -> 450 bp
  //   PA          (~2151 bp) -> 425 bp
  //   HA          (~1778 bp) -> 350 bp
  //   NP          (~1565 bp) -> 300 bp
  //   NA          (~1413 bp) -> 275 bp
  //   MP          (~1027 bp) -> 200 bp
  //   NS          ( ~890 bp) -> 175 bp
  // Override globally with params.max_skip_size, or accept per-segment defaults.
  // null = use segment-aware defaults; 0 = disabled.
  def dvg_skip_defaults = [
    'PB2': 450, 'PB1': 450, 'PA': 425,
    'HA':  350, 'NP':  300,
    'NA':  275,
    'MP':  200, 'NS':  175
  ]
  def max_skip = (params.max_skip_size != null)
      ? params.max_skip_size as int
      : dvg_skip_defaults.getOrDefault(segmentUpper, 175)

  // When true, unmapped reads are retained in the raw/unfiltered BAM by setting
  // the samtools -F flag to 0 (no filtering). Has no effect on the filtered BAM,
  // which always excludes unmapped reads via FLAGF_FILTERED.
  boolean output_unmapped = false
  if (params.output_unmapped_reads instanceof Boolean) {
    output_unmapped = params.output_unmapped_reads
  } else if (params.output_unmapped_reads != null) {
    output_unmapped = params.output_unmapped_reads.toString().toLowerCase() in ['1','true','yes','y']
  }

  // Master switch: enable dvg_filter.py post-alignment filtering.
  // When false, no filtering beyond FLAGF_RAW is applied.
  // The raw_bam is only written when filter_bam=true (run_filter=true).
  boolean filter_bam = false
  if (params.filter_bam instanceof Boolean) {
    filter_bam = params.filter_bam
  } else if (params.filter_bam != null) {
    filter_bam = params.filter_bam.toString().toLowerCase() in ['1','true','yes','y']
  }

  // Determine whether the filter branch should run:
  // requires nanopore platform, filter_bam=true, and at least one active filter criterion.
  boolean run_filter = (params.platform == 'nanopore') &&
                       filter_bam &&
                       (min_len > 0 || max_skip > 0 || max_clip_frac > 0)

  """
  MAP_OPTION=$map_option
  FLAGF_RAW=$flagF_raw
  FLAGF_FILTERED=$flagF_filtered
  MIN_LEN=$min_len
  MAX_SKIP=$max_skip
  MAX_CLIP_FRAC=$max_clip_frac
  OUTPUT_UNMAPPED=$output_unmapped
  FILTER_BAM=$filter_bam
  RUN_FILTER=$run_filter
  SEGMENT=$segmentUpper

  # ── Build minimap2 options ────────────────────────────────────────────────
  # Secondary alignments are ALWAYS allowed at the minimap2 level so that the
  # raw_bam retains multi-locus DVG mappings for breakpoint analysis.
  # Secondary suppression for the filtered BAM is applied via FLAGF_FILTERED.
  MM2_OPTS="-ax \$MAP_OPTION -t${task.cpus}"

  # ── Resolve raw BAM flag ──────────────────────────────────────────────────
  # If output_unmapped_reads=true, disable all samtools -F filtering (flag 0)
  # so unmapped reads are retained. Otherwise use FLAGF_RAW (default 4).
  SAMTOOLS_RAW_F=\$FLAGF_RAW
  if [ "\$OUTPUT_UNMAPPED" = "true" ]; then
    SAMTOOLS_RAW_F=0
  fi

  if [ "\$RUN_FILTER" = "true" ]; then
    echo "DVG filtering ENABLED: segment=\$SEGMENT, platform=${params.platform}, MIN_LEN=\$MIN_LEN, MAX_SKIP=\$MAX_SKIP, MAX_CLIP_FRAC=\$MAX_CLIP_FRAC, OUTPUT_UNMAPPED=\$OUTPUT_UNMAPPED" >&2

    # ── Pass 1: full alignment -> raw_bam ─────────────────────────────────
    # SAMTOOLS_RAW_F: 0 if output_unmapped=true, otherwise FLAGF_RAW (default 4).
    minimap2 \$MM2_OPTS $ref_fasta $reads \\
      | samtools view -h -@${task.cpus} -F \$SAMTOOLS_RAW_F \\
      | samtools view -@${task.cpus} -b \\
      | samtools sort -@${task.cpus} -O BAM -o $raw_bam

    samtools index $raw_bam

    # ── Pass 2: filtered BAM via dvg_filter.py ────────────────────────────
    # Reads from the already-written raw_bam to avoid the tee race condition.
    # FLAGF_FILTERED (default 2308) always excludes unmapped reads regardless
    # of output_unmapped_reads — unmapped reads have no role in variant calling.
    samtools view -h -@${task.cpus} -F \$FLAGF_FILTERED $raw_bam \\
      | dvg_filter.py \\
          --min-len \$MIN_LEN \\
          --max-skip \$MAX_SKIP \\
          --max-clip-frac \$MAX_CLIP_FRAC \\
          --segment \$SEGMENT \\
          --sample ${sample} \\
          --stats-out $dvg_stats \\
      | samtools view -@${task.cpus} -b \\
      | samtools sort -@${task.cpus} -O BAM > $bam

  else
    echo "DVG filtering DISABLED: segment=\$SEGMENT, platform=${params.platform}, MIN_LEN=\$MIN_LEN, MAX_SKIP=\$MAX_SKIP, MAX_CLIP_FRAC=\$MAX_CLIP_FRAC, OUTPUT_UNMAPPED=\$OUTPUT_UNMAPPED" >&2

    # ── No filter path: single pass directly to bam ──────────────────────
    # SAMTOOLS_RAW_F: 0 if output_unmapped=true, otherwise FLAGF_RAW (default 4).
    minimap2 \$MM2_OPTS $ref_fasta $reads \\
      | samtools view -h -@${task.cpus} -F \$SAMTOOLS_RAW_F \\
      | samtools view -@${task.cpus} -b \\
      | samtools sort -@${task.cpus} -O BAM > $bam
  fi

  samtools index $bam
  samtools stats  $bam > $stats
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
