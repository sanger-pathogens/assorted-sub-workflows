#!/usr/bin/env nextflow
include { BASECALLING } from './subworkflows/basecalling.nf'
def logo = NextflowTool.logo(workflow, params.monochrome_logs)

log.info logo

def printHelp() {
    NextflowTool.help_message("${workflow.ProjectDir}/schema.json", [],
    params.monochrome_logs, log)
}

def combine_metadata_maps(barcode, meta1, meta2, reads) {
    [meta1 + meta2, reads]
}

workflow {
    if (params.help) {
        printHelp()
        exit 0
    }

    Channel.fromPath(params.reference)
        .set{ reference }

    Channel.fromPath(params.target_regions_bed)
        .set{ target_regions_bed }

    Channel.fromPath(params.additional_metadata)
        .ifEmpty {exit 1, "${params.additional_metadata} appears to be an empty file!"}
        .splitCsv(header:true, sep:',')
        .map { meta -> ["${meta.barcode_kit}_${meta.barcode}", meta] }
        .set { additional_metadata }

    if (params.basecall) {
        raw_reads = Channel.fromPath("${params.raw_read_dir}/*.{fast5,pod5}", checkIfExists: true)
        BASECALLING(
            raw_reads,
            additional_metadata
        )
    }
  
workflow.onComplete {
    NextflowTool.summary(workflow, params, log)
}