#!/usr/bin/env nextflow
include { TRIMMING          } from './subworkflows/trimming.nf'
include { TR_FILTERING      } from './subworkflows/tr_filtering.nf'
include { HOST_READ_REMOVAL } from './subworkflows/host_read_removal.nf'
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

    if (params.run_trimmomatic){
        TRIMMING(decompressed_reads_ch)
        | set{ preprocessed_ch_1 }
    }
    else{
        preprocessed_ch_1 = decompressed_reads_ch
    }
    if (params.run_trf){
        TR_FILTERING(preprocessed_ch_1)
        | set{ preprocessed_ch_2 }
    }
    else{
        preprocessed_ch_2 = preprocessed_ch_1
    }
    if (params.run_bmtagger){
        HOST_READ_REMOVAL(preprocessed_ch_2)

        HOST_READ_REMOVAL.out.host_read_removal_out_ch
        | set{ preprocessed_ch_3 }
    }
    else{
        preprocessed_ch_3 = preprocessed_ch_2
    }

    COMPRESS_READS(preprocessed_ch_3)
    

    emit:
    preprocessed_reads_ch = COMPRESS_READS.out.compressed_reads_ch
    collated_host_reads_stats_ch = HOST_READ_REMOVAL.collated_host_reads_stats_ch
    collated_trimming_stats_ch = TRIMMING.collated_trimming_stats_ch


}