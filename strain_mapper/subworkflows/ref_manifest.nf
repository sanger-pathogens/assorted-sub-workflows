#!/usr/bin/env nextflow
//
// Check input samplesheet and get read channels
//
workflow REF_MANIFEST_PARSE {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    Channel
        .fromPath( samplesheet )
        .ifEmpty {exit 1, "Cannot find path file ${samplesheet}"}
        .splitCsv ( header:true, sep:',' )
        .map { create_ref_channels(it) }
        .filter{ meta, reference -> reference != null }  
        .set { references }

    emit:
    references // channel: [ val(meta), Path(reference) ]
}

// Function to get list of [ meta, fastq_1, fastq_2 ]
def create_ref_channels(LinkedHashMap row) {
    def meta = [:]
    meta.ID = row.ID

    Path reference = null

    def errors = []
    def array = []
    // check short reads
    if ( !(row.reference == 'NA') ) {
        if ( !file(row.reference).exists() ) {
            errors << "Reference fasta file does not exist:\n${row.reference}"
        }
        reference = file(row.reference)
    }

    if (errors) {
        errors.add(0, "Errors while parsing input manifest!")
        throw new Exception(errors.join('\n'))
    }

    array = [ meta, reference ]
    return array
}
