#!/usr/bin/env nextflow
include { TRIMMOMATIC       } from './modules/trimmomatic.nf'
include { TR_FILTERING      } from './subworkflows/tr_filtering.nf'
include { READ_REMOVAL      } from './subworkflows/read_removal.nf'
include { COMPRESS_READS       
          DECOMPRESS_READS  } from './modules/helper_processes.nf'

workflow PREPROCESSING {

    /*
    -----------------------------------------------------------------
    Preprocessing fastq files
    -----------------------------------------------------------------

    */

    take:
    reads_ch

    main:

    DECOMPRESS_READS(reads_ch)
    | set{ decompressed_reads_ch }

    if (params.run_trimmomatic && params.run_trf && params.run_bmtagger) {
        TRIMMOMATIC(decompressed_reads_ch)
        | TR_FILTERING
        | READ_REMOVAL
        | set{ processed_reads }
    } else if (params.run_trimmomatic && params.run_trf) {
        TRIMMOMATIC(decompressed_reads_ch)
        | TR_FILTERING
        | set{ processed_reads }
    } else if (params.run_trimmomatic && params.run_bmtagger) {
        TRIMMOMATIC(decompressed_reads_ch)
        | READ_REMOVAL
        | set{ processed_reads }
    } else if (params.run_trf && params.run_bmtagger) {
        TR_FILTERING(decompressed_reads_ch)
        | READ_REMOVAL
        | set{ processed_reads }
    } else if (params.run_trimmomatic) {
        TRIMMOMATIC(decompressed_reads_ch)
        | set{ processed_reads }
    } else if (params.run_trf) {
        TR_FILTERING(decompressed_reads_ch)
        | set{ processed_reads }
    } else if (params.run_bmtagger) {
        READ_REMOVAL(decompressed_reads_ch)
        | set{ processed_reads }
    } else {
        processed_reads = decompressed_reads_ch
    }

    COMPRESS_READS(processed_reads) 
    | set{ compressed_reads_ch }

    emit:
    compressed_reads_ch
}


def parse_mnf_meta(preprocessing_mnf) {
    // consensus_mnf <Channel.fromPath()>
    println(preprocessing_mnf)
    def mnf_ch =  Channel.fromPath(preprocessing_mnf)
                        | splitCsv(header: true, sep: ',')
                        | map {row -> 
                            // set meta
                            meta = [
                                id: row.sample_id,
                                sample_id: row.sample_id
                                ]

                            // set files
                            reads = [row.reads_1, row.reads_2]

                            // declare channel shape
                            tuple(meta, reads)
                        }
    return mnf_ch // tuple(index, [fastq_pairs])
}

// run preprocessing in isolation
workflow {
    mnf_ch = parse_mnf_meta(params.preprocessing_mnf)
    PREPROCESSING(mnf_ch)
}

def check_preprocessing_params(){
    /*
    -----------------------------------------------------------------
    Checks for necessary parameters and validates paths to ensure 
    they exist. Logs errors if any required parameters are missing.
    -----------------------------------------------------------------

    - **Output**: Number of errors encountered during the checks.

    -----------------------------------------------------------------

    */
    def errors = 0
    // was the kraken database provided?
    if (params.adapter_fasta == null){
        log.error("No adapter_fasta path provided")
        errors +=1
    }

    // if yes, is it a file which exists?
    if (params.adapter_fasta){
        adapter_fasta = file(params.adapter_fasta)
        if (!adapter_fasta.exists()){
            log.error("The adapter_fasta provided (${params.adapter_fasta}) does not exist.")
            errors += 1
        }
    }

    // check switchs

    if ((params.run_trimmomatic == false) && (params.run_trf == false) && (params.run_bmtagger == false)){
        log.error("All PREPROCESSING process switchs are off (run_trimmomatic = ${params.run_trimmomatic}; run_trf = ${params.run_trf} ; run_bmtagger = ${params.run_bmtagger}).")
    }

    return errors
}