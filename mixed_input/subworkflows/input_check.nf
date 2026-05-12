//
// Check input samplesheet and get read channels
//
include { FILTER_EXISTING_OUTPUTS       } from '../../pipeline_chaining/subworkflows/skip_downloaded.nf'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    Channel
        .fromPath( samplesheet )
        .ifEmpty {exit 1, "Cannot find path file ${samplesheet}"}
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channels(it) }
        .filter{ meta, reads_1, reads_2 -> reads_1 != null && reads_2 != null }  // Single end not supported
        .set { shortreads }
    
        if (params.only_new_input){
            FILTER_EXISTING_OUTPUTS(shortreads)
            FILTER_EXISTING_OUTPUTS.out.do_not_exist
            .set { shortreads_to_process }
        } else{
            shortreads_to_process = shortreads
        }

    emit:
    shortreads_to_process // channel: [ val(meta), file(reads_1), file(reads_2) ]
}

// Function to get list of [ meta, fastq_1, fastq_2 ]
def create_fastq_channels(LinkedHashMap row) {
    def meta = [:]
    meta.ID = row.ID

    Path fastq_1 = null
    Path fastq_2 = null

    def errors = []
    def array = []
    // check short reads
    if ( !(row.R1 == 'NA') ) {
        if ( !file(row.R1).exists() ) {
            errors << "R1 fastq file does not exist:\n${row.R1}"
        }
        fastq_1 = file(row.R1)
    }
    if ( !(row.R2 == 'NA') ) {
        if ( !file(row.R2).exists() ) {
            errors << "R2 fastq file does not exist:\n${row.R2}"
        }
        fastq_2 = file(row.R2)
    }
    
    if (errors) {
        errors.add(0, "Errors while parsing input manifest!")
        throw new Exception(errors.join('\n'))
    }

    array = [ meta, fastq_1, fastq_2 ]
    return array
}
