#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ILLUMINA } from './workflows/illumina'
include { NANOPORE } from './workflows/nanopore'
include { ASSEMBLIES } from './workflows/assemblies'

workflow NF_FLU {
    if (params.platform == 'illumina'){
        ILLUMINA ()
    } else if (params.platform == 'nanopore') {
        NANOPORE ()
    } else if (params.platform == 'assemblies') {
        ASSEMBLIES ()
    } else {
        exit 1, "ERROR: Unrecognized platform '${params.platform}'. Please choose illumina, nanopore, or assemblies."
    }
}

workflow {
    def json_schema = "$projectDir/nextflow_schema.json"

    if (params.help){
      def command = "nextflow run CFIA-NCFAD/nf-flu --input samplesheet.csv --platfrom <illumina/nanopore> samplesheet.csv -profile <singularity/docker/conda>"
      log.info NfcoreSchema.params_help(workflow, params, json_schema, command)
      exit 0
    }

    if (workflow.profile == 'slurm' && params.slurm_queue == "") {
      exit 1, "You must specify a valid SLURM queue (e.g. '--slurm_queue <queue name>' (see `\$ sinfo` output for available queues)) to run this workflow with the 'slurm' profile!"
    }

    if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
      custom_runName = workflow.runName
    }
    def summary_params = NfcoreSchema.params_summary_map(workflow, params, json_schema)
    log.info NfcoreSchema.params_summary_log(workflow, params, json_schema)

    NF_FLU ()

    workflow.onComplete = {
      // Log colors ANSI codes
      def c_reset = params.monochrome_logs ? '' : "\033[0m";
      def c_bold = params.monochrome_logs ? '' : "\033[1m";
      def c_red = params.monochrome_logs ? '' : "\033[0;31m";
      def c_green = params.monochrome_logs ? '' : "\033[0;32m";
      println """
      Pipeline execution summary
      ---------------------------
      Completed at : ${workflow.complete}
      Duration     : ${workflow.duration}
      Success      : ${c_bold}${workflow.success ? c_green : c_red}${workflow.success}${c_reset}
      Results Dir  : ${file(params.outdir)}
      Work Dir     : ${workflow.workDir}
      Exit status  : ${workflow.exitStatus}
      Error report : ${workflow.errorReport ?: '-'}
      """.stripIndent()
    }
    workflow.onError = {
        println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
    }
}

