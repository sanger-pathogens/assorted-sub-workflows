#!/usr/bin/env nextflow
include { TRIMMING       } from './subworkflows/trimming.nf'
include { TR_FILTERING      } from './subworkflows/tr_filtering.nf'
include { HOST_READ_REMOVAL      } from './subworkflows/host_read_removal.nf'
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
        TRIMMING(decompressed_reads_ch)
        | TR_FILTERING
        | HOST_READ_REMOVAL
        | set{ processed_reads }
    } else if (params.run_trimmomatic && params.run_trf) {
        TRIMMING(decompressed_reads_ch)
        | TR_FILTERING
        | set{ processed_reads }
    } else if (params.run_trimmomatic && params.run_bmtagger) {
        TRIMMING(decompressed_reads_ch)
        | HOST_READ_REMOVAL
        | set{ processed_reads }
    } else if (params.run_trf && params.run_bmtagger) {
        TR_FILTERING(decompressed_reads_ch)
        | HOST_READ_REMOVAL
        | set{ processed_reads }
    } else if (params.run_trimmomatic) {
        TRIMMING(decompressed_reads_ch)
        | set{ processed_reads }
    } else if (params.run_trf) {
        TR_FILTERING(decompressed_reads_ch)
        | set{ processed_reads }
    } else if (params.run_bmtagger) {
        HOST_READ_REMOVAL(decompressed_reads_ch)
        | set{ processed_reads }
    } else {
        processed_reads = decompressed_reads_ch
    }

    COMPRESS_READS(processed_reads) 
    | set{ compressed_reads_ch }
    

    emit:
    compressed_reads_ch
}
