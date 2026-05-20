# CFIA-NCFAD/nf-flu

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [[IRVC v1.2.0](https://github.com/PHAC-IRVC-Genomics/nf-flu/releases/tag/v1.2.0)] - 2026-05-XX

### Major Updates:

Re-implimented the Defective Viral Genome (DVG) filtering process as a standalone Python script (`dvg_filter.py`) called in `minimap2.nf`:

* Optional DVG filtering for Nanopore data (`--filter_bam`, default: `false`). When enabled, a two-pass alignment strategy produces both an unfiltered archive BAM and a filtered BAM suitable for variant calling. Filtering is handled by `dvg_filter.py` with four sequential CIGAR-based criteria and segment-aware thresholds:
  
     * Minimum aligned length (`--min_map_len`, default `null`: see below): minimum number of aligned (M, X and =) bases required to retain a read in the filtered BAM
     * Maximum single N/D operation (`--max_skip_size`, default `null`: see below): discards reads with a large internal deletion
     * Single contiguous alignment block: failsafe that discards reads with multiple alignment blocks not explained by a deletion (i.e., malformed/unusual CIGAR strings)
     * Soft-clip fraction (`--max_clip_frac`, default `0.15`): discards internally-mapping DVG reads with heavily clipped flanks larger than a typical primer overhang
       
* Both `--min_map_len` and `--max_skip_size` resolve to per-segment defaults at the Nextflow level in `minimap2.nf`:
     * `--min_map_len`, default: `null` = ~50% of segment length and `--max_skip_size`, default: `null` = ~20% segment length
          * Note: the `--max_skip_size` default (~20% segment length) should be comfortably larger than known biologically-important deletions (e.g., NA stalk deletions)
     * Both criteria can be overriden globally with a single value applied to every segment: (e.g., `--min_map_len 800` and `--max_skip_size 300`)
     * Both criteria can be disabled: (e.g., `--min_map_len 0` and `--max_skip_size 0`)
       
* Per-sample, per-segment filtering statistics and notable read identification (i.e., those with deletions in the biologically ambiguous range: 50 bp to `max_skip_size`) are written to dvg_filter_stats.json
  
* `--filter_bam`, `--min_map_len`, `--max_skip_size` and `--max_clip_frac` are configurable in `nextflow.config` and also from command line

* Commits [f44d861](https://github.com/PHAC-IRVC-Genomics/nf-flu/commit/f44d8613330a4966e414d67fa4d1b3729743445e), [af7a184](https://github.com/PHAC-IRVC-Genomics/nf-flu/commit/af7a184c84b72ea2515460c49320f88ca35e080c), [657e383](https://github.com/PHAC-IRVC-Genomics/nf-flu/commit/657e3834b8a01fc2d82bdb630da724543c21adf2), [fc09fde](https://github.com/PHAC-IRVC-Genomics/nf-flu/commit/fc09fde56f5c907473b47432afdf9f395a4e0fad), [593a845](https://github.com/PHAC-IRVC-Genomics/nf-flu/commit/593a8458fe3077fc6a2ad2120c02c9f9fb8e09d3), [59b4255](https://github.com/PHAC-IRVC-Genomics/nf-flu/commit/59b4255e454e467e85ac080a116dba38129d83bf), [28f25a8](https://github.com/PHAC-IRVC-Genomics/nf-flu/commit/28f25a89c65a1a30c8c0663f588bbcc184d0405d) and [1ec573b](https://github.com/PHAC-IRVC-Genomics/nf-flu/commit/1ec573bcc86461339f19b54e71c0f3bdf721bb98)

### Minor Updates:

* Update `Genin2` to v2.1.6 [[commit 547582b](https://github.com/PHAC-IRVC-Genomics/nf-flu/commit/547582b470575e96dec9902986c8fc3be2410bb3)]
  
* Update `Nextclade` to v3.21.0 [[commit 3621430 ](https://github.com/PHAC-IRVC-Genomics/nf-flu/commit/36214306bb43f07068e0b87603d08696b83cd14f)]
  
* Incorporated `cdsCoverage` from `Nextclade` results into `MultiQC` Report [[commit 3d86fa6](https://github.com/PHAC-IRVC-Genomics/nf-flu/commit/3d86fa68143fa9b923bde3eace0c9a4e9b5a71cf)]
  
* Updated `VADR` alert parameters for `assemblies` mode to pass the following alerts: `--alt_pass lowsim5s,lowsim3s,indf5pst,indf3pst` [[commit 1cafd1d](https://github.com/PHAC-IRVC-Genomics/nf-flu/commit/1cafd1dbc4eee4701ddc250b71433b54fcd92a2d)]
  
* Minor update to VADR output files to identify edgecase where segments pass VADR but could have short truncations [[commit e588c69](https://github.com/PHAC-IRVC-Genomics/nf-flu/commit/e588c69c853bbaa0b4ced0cedbe7ebe03006452d)]:
  * `vadr-annotation-alerts.txt` --> lists segments that fail VADR annotation
  * `vadr-annotation-alerts.txt` --> lists segments that pass VADR annotation but may have non-fatal issues such as 5' and/or 3' truncations
    
* Added error handling in `illumina.nf` and `nanopore.nf` workflows for read count calculation to skip and warn against corrupted .fastq files instead of failing the pipeline [[commit 1af430c](https://github.com/PHAC-IRVC-Genomics/nf-flu/commit/1af430c049af1c36bae9d2c892d1d9ce5fd15b88)] and [[commit 327b6f8](https://github.com/PHAC-IRVC-Genomics/nf-flu/commit/327b6f8cebfda0e569ad3720fede170ad0e67f15)]
  
* Updated `seqtk_seq.nf` to convert a degenerate nucleotide to one of its representative non-degenerate versions (similar to `Clair3`) in the selected reference sequence as `BCFTools` does not handle degenerate nucleotides during consensus generation, and will place an `N` instead of the appropriate nucleotide into the consensus sequence [[commit 0cfcc97](https://github.com/PHAC-IRVC-Genomics/nf-flu/commit/0cfcc97c37b242df2cdbb99af4e112bb6ca7149d)]

* Updated `subtyping_report.py` to backfill missing expected columns with empty strings in subtype results DataFrame to prevent KeyError causing the subtyping report process to fail [[commit 0d8ff12](https://github.com/PHAC-IRVC-Genomics/nf-flu/commit/0d8ff1235e09e2184c63afde0c6107cc169ec8fd)]

## [[IRVC v1.0.2](https://github.com/PHAC-IRVC-Genomics/nf-flu/releases/tag/v1.0.2)] - 2026-03-18
* Update IRMA to v1.2.0
* Updated Genoflu to v1.07
* Incorporated salient results from Nextclade analysis into MultiQC Report
* Minor code cleanup: removed legacy process `VCF_FILTER_FRAMESHIFT` from `modules_nanopore.config`; see [[CFIA-NCFAD/nf-flu PR-129](https://github.com/CFIA-NCFAD/nf-flu/pull/129/changes/cc2c1c9af7d94c9c02c178e2e8802e632d20313f)]

## [[IRVC v1.0.1](https://github.com/PHAC-IRVC-Genomics/nf-flu/releases/tag/v1.0.1)] - 2026-02-16
* Improved AWK-based filtering and removal of reads associated with defective viral genomes
    * Added criteria to reject reads with large internal deletions/skips and multi-block split alignments (i.e., short reads with secondary alignments)
    * New configurable parameters:
      * `max_gap_size`; pipeline default: 500 (bp); [corresponding minimap2 parameter: -G]
      * `allow_secondary`; pipeline default: false; [corresponding minimap2 parameter: --secondary=yes|no]
      * `max_skip_size`; pipeline default: 200 (bp); [Reads with N or D CIGAR operations larger than this threshold are rejected]
* Incorporated minor updates from CFIA-NCFAD v3.10.3
* Other minor QOL changes

## [[IRVC v1.0.0](https://github.com/PHAC-IRVC-Genomics/nf-flu/releases/tag/v1.0.0)] - 2025-12-11
* Forked from https://github.com/CFIA-NCFAD/nf-flu
* Modified default settings to IRVC standard
* Added AWK-based mapped read length BAM filtering for Nanopore sequencing
* Updated Genin2 to v2.1.5

## [[3.10.3](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.10.3)] - 2026-01-26

Update URLs (<https://api.figshare.com/v2/file/download/53449877> and <https://api.figshare.com/v2/file/download/53449874>) for test data used in GitHub Actions CI and NCBI Influenza sequences `ncbi_influenza_fasta` and metadata `ncbi_influenza_metadata` used by default by nf-flu.

## [[3.10.2](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.10.2)] - 2025-30-12

Update URL for test data used in GitHub Actions CI and NCBI Influenza sequences `ncbi_influenza_fasta` and metadata `ncbi_influenza_metadata` used by default by nf-flu.

## [[3.10.1](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.10.1)] - 2025-09-15

This patch release updates `nextflow_schema.json`, removing a duplicate `FLU` option from the `irma_module` ENUM. This fix is necessary in order to allow the pipeline to be launched using Seqera Platform.

## [[3.10.0](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.10.0)] - 2025-07-30

This minor release adds Genin2, updates Clair3 and fixes Nextclade consolidated output to use the correct Nextclade dataset tag/version.

### Changes

* feat: add Genin2 for European clade 2.3.4.4b H5Nx genotype prediction to the `illumina`, `nanopore` and `assemblies` workflows (i.e. `--platform illumina/nanopore/assemblies`).
* update: update Clair3 to 1.1.2 to fix a potential issue when using the `--enable_variant_calling_at_sequence_head_and_tail` option (#116)
* fix: Nextclade dataset tag/version used in consolidated tabular output file.
* dev: add `run-assemblies-test.sh` and update `run-*-test.sh` scripts in `tests/`. Better default `--outdir` for each test script and downloading of VADR model prior to running workflow to allow sharing of model tar.gz between different tests.

## [[3.9.0](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.9.0)] - 2025-04-07

This minor release adds Nextclade analysis of assembled Influenza genome sequences against 30 Nextclade Influenza-related datasets by default and updates the Influenza sequences used by nf-flu (downloaded from NCBI 2025-04-04; 809,739 unique sequences and metadata).

The specific Nextclade datasets and optionally versions (tags) can be configured with a headerless CSV file. Nextclade results are aggregrated across samples and datasets and filtered for positive results into a single Nextclade TSV (tab-separated values) report with additional fields capturing sample, dataset name and dataset version/tag information as well as Nextclade and Nextclade dataset specific results.

### Changes

* update: Influenza sequences and metadata from NCBI (2025-04-04). 809,739 non-redundant, unique sequences were retrieved along with their metadata. Added documentation for how to update Influenza sequences for use with nf-flu (see [docs/update_seqs_db.md](docs/update_seqs_db.md))
* feat: added Nextclade (v3.12.0) analysis subworkflow against 30 Influenza-related Nextclade datasets with a convenient aggregation and summarization of useful results into a single Nextclade TSV report.
* update: GenoFLU 1.05 -> 1.06 (#112)

## [[3.8.1](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.8.1)] - 2025-03-26

This patch release updates Clair3, the variant caller for Nanopore sequence data, to 1.0.11. Clair 1.0.11 adds an option to enable variant calling at the head and tail 16bp of each sequence (`--enable_variant_calling_at_sequence_head_and_tail`). This option is enabled by default in the nf-flu workflow to ensure that the 16bp at the start and end of each of the 8 segments of IAV and IBV are variant called. It should be noted that the developers of Clair3 note that results are used "with caution because alignments are less reliable in the regions, and there would be insufficient context to be fed to the neural network for reliable calling".

A minor issue with the MultiQC report was also fixed where sample names were not cleaned properly. The `.bcftools_filt` extension was added to `extra_fn_clean_exts` in `assets/multiqc_config.yaml`.

### Changes

* fix: MultiQC report sample name cleaning. Added `.bcftools_filt` to `extra_fn_clean_exts` in `assets/multiqc_config.yaml`.
* update: Clair3 1.0.10 -> 1.0.11
* fix: Clair3 not variant calling the ends of each segment enable variant calling at the head and tail 16bp of each sequence (`--enable_variant_calling_at_sequence_head_and_tail`) (#61)
* dev: move Clair3 arguments and options to `conf/modules_nanopore.config`. This should allow users to change Clair3 options more easily using custom Nextflow config files (e.g. `nextflow run CFIA-NCFAD/nf-flu -c clair3-custom.config ...`).
* test: added nf-test for `clair3.nf` to with simulated test data for head and tail variant calling with the `--enable_variant_calling_at_sequence_head_and_tail` option.

## [[3.8.0](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.8.0)] - 2025-03-25

This release adds the `--platform assemblies` mode for analysis of FASTA sequences along with `--input /path/to/fasta-dir/` to specify the directory containing the FASTA sequences.

### Changes

* feat: analysis of previously assembled IAV FASTA sequences with the addition of a new analysis mode via `--platform assemblies`. Use along with `--input /path/to/fasta-dir/` to specify the directory containing the FASTA sequences.
* fix: `bin/cleavage_site.py` short cleavage site index access error (#106)
* fix: `cleavage_site.nf` version output issue (#105)
* fix: low abundance indels appearing in consensus sequences despite major/minor allele frequency thresholds. Explicitly excluding non-SNP variants below the major allele fraction prior to consensus sequence generation with Bcftools consensus.
* fix: subtyping report issue with some poor quality IBV sequences (#107)
* dev: add nf-test for VCF filtering and consensus sequence generation from VCF with low AF indels.
* dev: replaced `vcf_filter_frameshift.py` with Bcftools filter commands.

## [[3.7.0](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.7.0)] - 2025

This minor release adds GenoFLU for H5 genotyping and HA cleavage site output with VADR annotations. This release also adds a script to classify HA cleavage sites based on mono-/multibasicity and low/high pathogenicity.

### Changes

* feat: GenoFLU v1.05 for H5 genotyping.
* feat: Added `--custom_flu_minfo` option to specify custom `flu.minfo` for VADR. The default `flu.minfo` is the same as the VADR flu v1.6.3-2 model except that it includes cleavage site info. Feature table, GenBank and GFF files should now have a `misc_feature` for HA cleavage site info.
* feat: `bin/cleavage_site.py` to classify HA cleavage sites.
* feat: Added VADR subtype prediction into subtyping report. VADR subtype predictions are pulled from the output `.mdl` files.
* feat: Added subtyping report output directory containing CSV for each sheet in the Excel report.
* fix: MultiQC converts the general info table into a violin plot if there are more than 500 rows in the table by default. Added `max_table_rows: 1000000` to `multiqc_config.yaml` to avoid this conversion in most cases.

## [[3.6.2](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.6.2)] - 2025-01-07

This patch release fixes issues relating to subtype prediction (N5) (#100), Apptainer usage (#95) and IRMA read length threshold (#99).

### Changes

* chore: renamed: `bin/parse_influenza_blast_results.py` -> `bin/subtyping_report.py`
* fix: N5 sequences being typed as N1 due to the high proportion of lower % identity hits to N1 sequences (#100). In `bin/subtyping_report.py`, H and N subtype is predicted based on determining what the subtype is using the BLAST analysis results starting at a % identity threshold of 99% and decrementing by 1% until a subtype or the minimum % identity is reached (default: 85%). At least 3 hits are required to determine a subtype at a particular threshold. If no subtype is determined, the subtype is set to "N/A".
* fix: added back missing results columns to subtyping report H and N subtyping sheets.
* fix: Added workflow parameter `--irma_min_len` to be able to change the minimum read length threshold for IRMA assembly (`MIN_LEN`) and set default to 50 instead of 125. nf-flu should now be compatible with BGI sequencing data producing shorter paired-end reads by running with `--platform illumina` (#99).
* fix: `-profile apptainer` is functionally the same as `-profile singularity`. The same configuration is set for the Apptainer profile as for the Singularity one. If a user has Apptainer installed, running `$ singularity ...` and `$ apptainer ...` should be equivalent, e.g. both `$ apptainer --version` and `$ singularity --version` produce `apptainer version 1.3.6`. (#95)
* ci: updated ci.yml for better cache handling and inter-job caching of VADR flu model tar.gz
* config: `--max_top_blastn` default changed 3 -> 5. Top 5 BLASTN hits will be shown for each segment for each sample in subtyping report.

## [[3.6.1](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.6.1)] - 2024-12-13

This patch release fixes an issue with Clair3 not producing variant calls for some regions due to full-alignment not being triggered. This issue was resolved by adding `--var_pct_phasing=1`, `--var_pct_full=1` and `--ref_pct_full=1` to the Clair3 command line.

### Changes

* fix: Added `--var_pct_phasing=1`, `--var_pct_full=1` and `--ref_pct_full=1` to Clair3 command line to ensure full-alignment is triggered for all reads to avoid missing variant calls in some regions.
* fix: Added `stageAs: "input*/*"` to `CAT_NANOPORE_FASTQ` process input channels to ensure that input files are not concatenated with themselves in an infinite loop until disk space is exhausted in rare cases.
* feat: Don't save NCBI Influenza reference sequences, metadata CSV and BLAST DB to the output directory by default. Added `--save_ncbi_db` and `--save_blastdb` workflow params to save these files to the output directory if desired.
* docs: Updated README.md to mention Apptainer. Updated `usage.md` to describe new workflow params. Updated `output.md` to better describe BLAST subtyping results.

## [[3.6.0](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.6.0)] - 2024-12-02

This minor release adds [FluMut](https://github.com/izsvenezie-virology/FluMut) to "to search for molecular markers with potential impact on the biological characteristics of Influenza A viruses of the A(H5N1) subtype."

### Changes

* **feat**: Added FluMut (v0.6.3)

## [[3.5.3](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.5.3)] - 2024-11-01

This patch release fixes an issue ([#22](https://github.com/peterk87/nf-flu/issues/22)) with Illumina paired-end read analysis by IRMA producing empty consensus sequences when the forward and reverse reads do not contain "1:N:0:." or "2:N:0:." in the FASTQ header lines.

## [[3.5.2](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.5.2)] - 2024-10-18

This patch release fixes a few issues when running the pipeline.

### Changes

* fix: better handling of empty IRMA consensus sequences to avoid downstream analysis errors with VADR and BLASTN ([peterk87/nf-flu #22](https://github.com/peterk87/nf-flu/issues/22))
* fix: Clair3 `versions.yml` indentation issue ([#87](https://github.com/CFIA-NCFAD/nf-flu/issues/87))
* fix: removed capturing of cat and gzip versions in CAT_ILLUMINA_FASTQ process ([#46](https://github.com/CFIA-NCFAD/nf-flu/issues/46)) to avoid issue in some execution environments.
* docs: update README.md

## [[3.5.1](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.5.1)] - 2024-10-08

This patch release fixes an issue ([#84](https://github.com/CFIA-NCFAD/nf-flu/issues/84)) with long sample names (over 50 characters) causing VADR to fail. `--noseqnamemax` has been added to the default arguments for VADR to avoid this issue.

### Changes

* **fix**: Added `--noseqnamemax` to VADR default arguments to avoid issues with long sample names causing VADR to fail.
* **config**: Output directory paths for IRMA and Bcftools consensus VADR annotation results were made more explicit and clear for the Illumina workflow.

## [[3.5.0](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.5.0)] - 2024-10

This release expands the Illumina workflow by adding BLAST analysis, coverage plots, variant calling, and MultiQC reports. Modifications were made to existing modules, and new modules were added.

### Changes

* **feat**: Added variant calling, BLAST analysis, coverage plots, and MultiQC to the Illumina workflow to match the capabilities of the Nanopore workflow.
* **feat**: Introduced a new module, Freebayes, for Illumina variant calling.
* **refactor**: Rearranged the Illumina workflow to integrate the new changes and enhance compatibility.
* **update**: Updated Bcftools filtering to add missing tags with `fill-tags` plugin and to set genotype with the `setGT` plugin based on `major/minor_allele_fraction` thresholds to influence consensus sequence output.
* **config**: Changed process labels for IRMA and MultiQC modules to "long" to avoid timeouts for large short-read datasets.
* **enhance**: Changed VADR staged file to use the FTP NCBI link to bypass certificate issues during Nextflow staging.
* **rollback**: Reverted VADR containers to an earlier version to resolve potential issues on Singularity.
* **refactor**: Rearranged `modules_illumina.config` for consistency with the updated workflow.
* **container**: Switched to Biocontainers images for Clair3 v1.0.10. [Issue](https://github.com/HKU-BAL/Clair3/issues/98) with full alignment not working with the Biocontainers Docker/Apptainer images seems to have been resolved. This should also resolve an issue with CI where it would fail due to not being able to pull the official Clair3 image [hbukal/clair3](https://hub.docker.com/r/hkubal/clair3) from Docker Hub.
* **dev**: Added `tests/run-illumina-test.sh` to make it more convenient to run the Illumina test locally with the same conditions as GitHub Actions CI.

## [[3.4.1](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.4.1)] - 2024-08-02

This patch release fixes an issue (#75) with CAT_ILLUMINA_FASTQ where `1:N:0:.` or `2:N:0:.` may be mistakenly appended
to Q-score lines beginning with `@`.

### Changes

* fix: updated Perl regex to better match Illumina FASTQ header lines starting with `@`. At least one space ` ` is expected in the header line. Match regex has been changed to `/^@.* .*/` from `/^@.*/` so hopefully Q-score lines should not be matched anymore.
* dev: replaced nf-core/modules `DUMPSOFTWAREVERSIONS` with [mqc_versions_table v0.2.0](https://github.com/CFIA-NCFAD/nim-mqc-versions-yml/releases/tag/0.2.0) Nim statically compiled binary to parse `versions.yml` and output necessary YAML with HTML content for display of process and tool versions table in MultiQC report. In theory `DUMPSOFTWAREVERSIONS` should be using the same Docker/Singularity image/Conda env as the MultiQC process, but `DUMPSOFTWAREVERSIONS` uses an older version of MultiQC and only uses it for the pyyaml library. `mqc_versions_table` was developed to handle this instead with a small 200KB binary instead.
* dev: harmonize Docker/Singularity containers and Conda envs used across processes.
* ci: use `symlink` mode for `publishDir` by default for `test_nanopore.config` and `test_illumina.config` to limit disk usage during CI.
* Updated to Bioconda channel VADR v1.6.4 since the STAPH-B offered container with the flu model packaged is very large at 6GB vs 1.45GB for quay.io/biocontainers/vadr 1.6.4. However, it's now required that the flu model be downloaded and installed prior to VADR annotation with `--vadr_model_targz`. The default model tarball is the `vadr-models-flu-1.6.3-2.tar.gz` (38MB) from the NCBI FTP site uploaded to [Zenodo](https://zenodo.org/records/13261208).

## [[3.4.0](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.4.0)] - 2024-07-24

This release adds Influenza virus sequence annotation using VADR.

### Changes

* Add VADR for Influenza consensus sequence annotation
* Add table2asn for Feature Table conversion to Genbank
* Add pre- and post-table2asn processing to workaround sequence ID length limits imposed by table2asn when converting from Feature Table format to Genbank

## [[3.3.10](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.3.10)] - 2024-05-31

Fix MultiQC report generation due to module filter paths not working like in v1.12.

### Software Updates

* multiqc: `1.21` -> `1.22.1`

### Changes

* test: add `tests/run-nanopore-test.sh` to conveniently run Nanopore test locally

## [[3.3.9](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.3.9)] - 2024-05-30

Long overdue software updates release.

### Software Updates

* bcftools: `1.15.1` -> `1.20`
* blast: `2.14.0` -> `2.15.0`
* clair3: `1.0.5` -> `1.0.9`
* minimap2: `2.24` -> `2.28`
* mosdepth: `0.3.3` -> `0.3.8`
* multiqc: `1.12` -> `1.21`
* seqtk: `1.3` -> `1.4`

### Changes

* dev: update GitHub Actions versions for CI and linting workflows

## [[3.3.8](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.3.8)] - 2024-02-16

This bugfix patch release fixes an issue where a large number of ambiguous bases in the IRMA consensus can hinder
reference selection (#67). This release also addresses an issue with using the Clair3 Biocontainers image resulting in
incomplete variant calling results, affecting nf-flu executions with the `docker` or `singularity` profiles. The
official Clair3 image is used instead. nf-flu executions using Conda and Mamba are unaffected.

### Changes

* Create majority consensus from IRMA `allAlleles.txt` files for BLASTN search
* Add `irma-alleles2fasta.v`, statically compiled binary (`irma-alleles2fasta`) and Bash build script for parsing IRMA
  `allAlleles.txt` to output naive majority consensus (i.e. whatever the top non-dash allele is at each position) so
  that the sequence used for BLASTN search does not contain any ambiguous bases.
* Updated nanopore.nf subworkflow to use IRMA majority consensus with no ambiguous bases for BLASTN search so that
  longer more contiguous matches are possible to aid in top reference sequence selection in some cases.
* Updated parse_influenza_blast_results.py to better handle extraction of sample name and segment number from BLASTN
  query accession/version (qaccver).
* Using official Clair3 Docker image and updating Clair3 to v1.0.5

## [[3.3.7](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.3.7)] - 2024-02-09

This bugfix patch release fixes an issue with mislabeling of PB1 and PB2 segments for Influenza B virus results (#65).

## [[3.3.6](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.3.6)] - 2023-11-01

### Fixes

* docs updated to show proper profile to run test profiles for Illumina and Nanopore locally (#52)
* `test_nanopore` profile has been updated to run locally with [the test samplesheet.csv updated with URLs to FASTQ files at CFIA-NCFAD/nf-test-datasets](https://github.com/CFIA-NCFAD/nf-test-datasets/blob/nf-flu/samplesheet/samplesheet_test_nanopore_influenza.csv)
* read samplesheet CSV in `parse_influenza_blast_results.py` with all columns read as string rather than inferred (#54)
* handle cloud storage paths and non-HTTP/FTP URLs in user samplesheets (#55)

## [[3.3.5](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.3.5)] - 2023-09-15

### Fixes

* handling of empty IRMA `amended_consensus/` when running a negative control or blank sequence (#47)

## [[3.3.4](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.3.4)] - 2023-08-18

### Fixes

* Subtyping report summary sheet "1_Subtype Predictions" shows only N subtype results

## [[3.3.3](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.3.3)] - 2023-08-16

This release fixes issues with subtype report generation script (`parse_influenza_blast_results.py`), primarily subtype predictions being `N/A` for samples where the top BLAST hits are user-specified sequences for the HA and NA segments.

### Fixes

* subtype prediction based off majority H/N prediction of all BLAST hits instead of just the top X matches (#40)
* the top hit for H/N can also be a user-specified sequence without subtype information
* top segment matches are now sorted by sample name, segment name and BLAST bitscore
* output concatenated Nanopore FASTQ to `${outdir}/fastq` by default (#43)
* Handle ambiguous bases in reference sequences by having Clair3 not convert those positions to N and Bcftools produce a warning instead of an error (#42)

### Changes

* subtyping report results are now ordered in the same order as the input `samplesheet.csv`, that is the order of the samples in the report is the same as the order of the samples in the `samplesheet.csv` file

## [[3.3.2](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.3.2)] - 2023-08-03

This patch release fixes an IBV subtype/genotype parsing issue when generating subtyping report using the new metadata format introduced in 3.3.0 ([#32](https://github.com/CFIA-NCFAD/nf-flu/issues/32)).

## [[3.3.1](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.3.1)] - 2023-08-02

### Fixes

* Conda/Mamba env creation when using `conda`/`mamba` profile (#35)

## [[3.3.0](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.3.0)] - 2023-07-11

This release migrates to more recently updated Influenza virus sequences since the last update for the [NCBI Influenza DB FTP data](https://ftp.ncbi.nih.gov/genomes/INFLUENZA/) was in 2020-10-13. By default, all Orthomyxoviridae virus sequences were parsed from the daily updated NCBI Viruses [`AllNucleotide.fa`](https://ftp.ncbi.nlm.nih.gov/genomes/Viruses/AllNucleotide/) and [`AllNuclMetadata.csv.gz`](https://ftp.ncbi.nlm.nih.gov/genomes/Viruses/AllNuclMetadata/AllNuclMetadata.csv.gz) and uploaded to [Figshare](https://figshare.com/articles/dataset/2023-06-14_-_NCBI_Viruses_-_Orthomyxoviridae/23608782) as Zstd compressed files. nf-flu no longer uses the [influenza.fna.gz](https://ftp.ncbi.nih.gov/genomes/INFLUENZA/influenza.fna.gz) and [genomeset.dat.gz](https://ftp.ncbi.nih.gov/genomes/INFLUENZA/genomeset.dat.gz) files for Influenza sequences and metadata, respectively.

### Fixes

* More up-to-date Influenza sequences database used by default (#24)

## [[3.2.1](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.2.1)] - 2023-07-07

### Fixes

* Empty BLAST results file parsing `NoDataError` (#27) (Thanks @MatFish for reporting this issue!)

## [[3.2.0](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.2.0)] - 2023-06-22

### Added

* Influenza B virus support (#14)
* Polars for faster parsing of BLAST results (#14)

### Fixes

* Irregular Illumina paired-end FASTQ files not producing IRMA assemblies (#20)

### Updates

* Updated README.md to include references and citations

## [[3.1.6](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.1.6)] - 2023-05-31

This is a patch release for a minor change to use Biocontainers Docker and Singularity images for Clair3 to avoid hitting limits on pulls from Docker Hub and since Biocontainers images are half the size of [hkubal/clair3](https://hub.docker.com/r/hkubal/clair3/) images.

Also, updated CI workflow and added issue template forms for feature request and questions.

## [[3.1.5](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.1.5)] - 2023-05-30

### Added

* `--use_mamba` to enable using [Mamba](https://github.com/mamba-org/mamba/) in place of Conda when using `-profile conda` for faster creation of Conda environments

### Updates

* Clair3: 0.1.10 -> 1.0.2

### Fixes

* user-specified Clair3 models not being found ([#11](https://github.com/CFIA-NCFAD/nf-flu/issues/11))
* Conda profile not enabling Conda ([#15](https://github.com/CFIA-NCFAD/nf-flu/issues/15))
* IRMA wanting too much `/tmp` space; IRMA's tmp dir will be output to the current working directory of the process job ([#13](https://github.com/CFIA-NCFAD/nf-flu/issues/13)) (Thanks @Codes1985 for reporting and solving this issue!)

## [[3.1.4](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.1.4)] - 2023-05-17

This release addresses issue [#11](https://github.com/CFIA-NCFAD/nf-flu/issues/11) adding a new option `--clair3_user_variant_model <PATH TO CLAIR3 MODEL>` to allow user to provide a Clair3 model not included with Clair3, e.g. a [Rerio](https://github.com/nanoporetech/rerio) Clair3 model for r10 flowcells.

## [[3.1.3](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.1.3)] - 2023-04-28

Patch release to fix issue to handle lowercase subtypes (e.g. `h1n5`) from NCBI Influenza DB.

## [[3.1.2](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.1.2)] - 2022-09-01

Patch release to fix issue when user reference sequences FASTA specified, but Channel from file is not treated as a value. Code has been reverted to use `file` Nextflow function.

## [[3.1.1](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.1.1)] - 2022-08-31

Patch release to fix issue when a user-specified sequences FASTA is provided and the FASTA is concatenated with the NCBI influenza sequences FASTA, but there is no new-line character at the end of the FASTA files. New line characters are added to the FASTA files to avoid incorrect concatenation.

## [[3.1.0](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.1.0)] - 2022-05-31

The workflow's name has been changed from `nf-iav-illumina` to `nf-flu` and the official repo for `nf-flu` will be [CFIA-NCFAD/nf-flu](https://github.com/CFIA-NCFAD/nf-flu/) going forward.

* Added back `bin/fastq_dir_to_samplesheet.py` for Illumina `--input` samplesheet creation from Illumina FASTQ reads directory
* Fixed [issue #12](https://github.com/peterk87/nf-flu/issues/12). Nanopore sample sheet can specify a mix of single FASTQ files and/or directories containing FASTQ files. Different reads with the same sample name will be merged prior to analysis. FASTQs can be GZIP compressed and have the extensions: `.fastq`, `.fq`, `.fastq.gz`, `.fq.gz`. Updated CI tests to test for this flexible sample sheet handling.
* Switched to GitHub YAML form for bug report template from Markdown template.
* CI tests now output `results/pipeline_info/` and `.nextflow.log` as artifacts for easier debugging of issues.

## [[3.0.0](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.0.0)] - 2022-05-24

This is a major release adding a Nanopore influenza sequence analysis subworkflow using IRMA for initial assembly and BLAST against NCBI Influenza DB sequences and optionally, user-specified sequences to identify the top reference sequence for each segment for each sample. A standard read mapping/variant calling analysis is performed: for each sample, Nanopore reads are mapped separately against each gene segment reference sequence using Minimap2; variant calling of read alignments is performed using Clair3; depth-masked consensus sequence is generated using Bcftools. Consensus sequences are BLAST searched against NCBI Influenza (and user-specified sequences) to generate a BLAST summary report and H/N subtyping report. MultiQC is used to summarize results into an interactive HTML report.

NOTE: Read mapping/variant calling analysis has not been ported to the Illumina sequence analysis subworkflow.

## [[2.0.1](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/2.0.1)] - 2021-06-15

Patch release to fix issue [#5](https://github.com/CFIA-NCFAD/nf-flu/issues/5); added check that IRMA `amended_consensus/` exists before concatenation of consensus FASTA files.

## [[2.0.0](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/2.0.0)] - 2021-06-10

### :warning: Major enhancements

* Samplesheet input (`--input samplesheet.csv`) replaces path to reads (`--reads "reads/*_R{1,2}_*.fastq.gz"`). Sample sheet can be tab-delimited (TSV) or CSV and must have a header line and 3 columns (sample name, FASTQ path/URL to forward reads, FASTQ path/URL to reverse reads).
* Pipeline has been re-implemented in [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html)
* All software containers are now exclusively obtained from [Biocontainers](https://biocontainers.pro/#/registry)
* Updated minimum Nextflow version to `v21.04.0` (see [nextflow#572](https://github.com/nextflow-io/nextflow/issues/1964))
* Add IRMA params
  * `irma_module`: IRMA module (default: `FLU-utr`)
  * `keep_ref_deletions`: set consensus sequence deletion by ambiguation (i.e. replace ref seq with Ns) (default: `true`)
* Add BLAST subtyping params:
  * `pident_threshold`: % identity threshold (default: `0.85`)
  * `min_aln_length`: min alignment length (default: `50`)
* Replace Azure Pipelines CI with GitHub Actions CI
* add `nextflow_schema.json` and nf-core helper Jar file and Groovy scripts for params validation, printing help
* Use nf-core modules where possible
* Use nf-core module style for all processes
* Added usage and output docs
* Updated README

### Parameters

| Old parameter | New parameter                         |
|:--------------|:--------------------------------------|
| `--reads`     | `--input`                             |
|               | `--irma_module`                       |
|               | `--keep_ref_deletions`                |
|               | `--pident_threshold`                  |
|               | `--min_aln_length`                    |
|               | `--ncbi_influenza_fasta`              |
|               | `--ncbi_influenza_metadata`           |
|               | `--slurm_queue_size`                  |
|               | `--publish_dir_mode`                  |
|               | `--validate_params`                   |
|               | `--enable_conda`                      |
|               | `--singularity_pull_docker_container` |
|               | `--show_hidden_params`                |
|               | `--schema_ignore_params`              |

* **NB:** Parameter has been **updated** if both old and new parameter information is present.
* **NB:** Parameter has been **added** if just the new parameter information is present.
* **NB:** Parameter has been **removed** if new parameter information isn't present.

### Software dependencies

Note, since the pipeline is now using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Dependency | Old version | New version |
|:-----------|:------------|:------------|
| `blast`    | 2.9.0       | 2.10.0      |
| `irma`     | 0.6.7       | 1.2.1       |
| `python`   | 3.7.3       | 3.9.0       |

* **NB:** Dependency has been **updated** if both old and new version information is present.
* **NB:** Dependency has been **added** if just the new version information is present.
* **NB:** Dependency has been **removed** if new version information isn't present.
