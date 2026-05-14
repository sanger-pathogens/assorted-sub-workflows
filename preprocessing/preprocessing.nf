#!/usr/bin/env nextflow
include { TRIMMING          } from './subworkflows/trimming.nf'
include { TR_FILTERING      } from './subworkflows/tr_filtering.nf'
include { HOST_READ_REMOVAL } from './subworkflows/host_read_removal.nf'
include { COMPRESS_READS;
          COMPRESS_SINGLE;       
          DECOMPRESS_READS  } from './modules/helper_processes.nf'

workflow PREPROCESSING {

    /*
    -----------------------------------------------------------------
    Preprocessing fastq files
    -----------------------------------------------------------------

    */

    take:
    reads_ch //meta, read_1, read_2

    main:

    DECOMPRESS_READS(reads_ch)
    | set{ decompressed_reads_ch }

    if (params.run_trimmomatic){
        TRIMMING(decompressed_reads_ch)

        preprocessed_ch_1          = TRIMMING.out.paired_ch
        collated_trimming_stats_ch = TRIMMING.out.stats_ch
        unpaired_ch_1              = TRIMMING.out.unpaired_ch
    }
    else {
        preprocessed_ch_1          = decompressed_reads_ch
        collated_trimming_stats_ch = Channel.empty()
        unpaired_ch_1              = Channel.empty()
    }

    if (params.run_trf){
        TR_FILTERING(preprocessed_ch_1, unpaired_ch_1)

        preprocessed_ch_2   = TR_FILTERING.out.trf_paired_ch
        unpaired_ch_2       = TR_FILTERING.out.trf_unpaired_ch
    }
    else {
        preprocessed_ch_2   = preprocessed_ch_1
        unpaired_ch_2       = unpaired_ch_1
    }

    if (params.run_bmtagger){
        HOST_READ_REMOVAL(preprocessed_ch_2, unpaired_ch_2)

        preprocessed_ch_3            = HOST_READ_REMOVAL.out.paired_ch
        unpaired_ch_3                = HOST_READ_REMOVAL.out.unpaired_ch
        collated_host_reads_stats_ch = HOST_READ_REMOVAL.out.collated_host_reads_stats_ch
    }
    else {
        preprocessed_ch_3            = preprocessed_ch_2
        unpaired_ch_3                = unpaired_ch_2
        collated_host_reads_stats_ch = Channel.empty()
    }

    COMPRESS_READS(preprocessed_ch_3)

    COMPRESS_SINGLE(unpaired_ch_3)

    emit:
    preprocessed_reads_ch = COMPRESS_READS.out.compressed_reads_ch
    collated_trimming_stats = collated_trimming_stats_ch
    collated_host_reads_stats = collated_host_reads_stats_ch
    unpaired_reads_compressed_ch = COMPRESS_SINGLE.out.compressed_reads_ch
}