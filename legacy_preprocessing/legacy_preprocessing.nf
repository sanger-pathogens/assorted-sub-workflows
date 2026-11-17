#!/usr/bin/env nextflow
include { TRIMMING          } from './subworkflows/trimming.nf'
include { TR_FILTERING      } from './subworkflows/tr_filtering.nf'
include { HOST_READ_REMOVAL } from './subworkflows/host_read_removal.nf'
include { COMPRESS_READS       
          DECOMPRESS_READS  } from './modules/helper_processes.nf'

workflow PREPROCESSING_ILLUMINA {

    /*
    -----------------------------------------------------------------
    Preprocessing fastq files
    -----------------------------------------------------------------

    */

    take:
    reads_ch

    main:

    // Illumina preprocessing subworkflow, passing only Illumina reads into the preprocessing steps
    if (params.preprocessing) {
        DECOMPRESS_READS(reads_ch)
        | set{ decompressed_reads_ch }

        if (params.run_trimmomatic){
            TRIMMING(decompressed_reads_ch)

            TRIMMING.out.trimmed_fastqs
            | set{ preprocessed_ch_1 }

            TRIMMING.out.collated_trimming_stats_ch
            | set { collated_trimming_stats_ch }
        }
        else{
            preprocessed_ch_1 = decompressed_reads_ch
            collated_trimming_stats_ch = Channel.empty()
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
            
            HOST_READ_REMOVAL.out.collated_host_reads_stats_ch
            | set { collated_host_reads_stats_ch }
        }
        else{
            preprocessed_ch_3 = preprocessed_ch_2
            collated_host_reads_stats_ch = Channel.empty()
        }

        COMPRESS_READS(preprocessed_ch_3)
        | set{ preprocessed_reads_ch }

    } else {
        // for when turning off preprocessing or for non-illumina data, 
        // just pass through the reads without preprocessing and empty channels for stats
        // so to streamline downstream workflow compatibility
        preprocessed_reads_ch = reads_ch
        collated_trimming_stats_ch = Channel.empty()
        collated_host_reads_stats_ch = Channel.empty()
    }
    emit:
    preprocessed_reads_ch
    collated_trimming_stats_ch
    collated_host_reads_stats_ch
}

workflow PREPROCESSING {
    /*
    -----------------------------------------------------------------
    Preprocessing fastq files
    -----------------------------------------------------------------

    */

    take:
    reads_ch

    main:
    // this gating might be redundant given that IRODS_EXTRACTOR already outputs platform-segregated channels, 
    // but it is a good safety check to ensure that only Illumina reads are passed into the preprocessing steps
    // items are [meta, read_1, read_2] (paired) or [meta, reads] (ONT), so only meta is read here
    reads_ch.branch{ reads ->
        // PUBLIC_PAIRED: paired reads from ENA (see MIXED_INPUT), preprocessed like Illumina
        illumina: reads[0].Platform in ["ILLUMINA", "PUBLIC_PAIRED"]

        ONT: reads[0].Platform == "ONT"

        other: true
    }
    | set { reads_ch }

    // only Illumina reads are preprocessed for now so just make a not processed channel
    reads_ch.ONT
    | mix(reads_ch.other)
    | set { not_preprocessed_ch }

    not_preprocessed_ch.subscribe { reads ->
        log.info "${reads[0].ID}: passed through without preprocessing, as only Illumina reads are preprocessed (Platform: ${reads[0].Platform})"
    }

    PREPROCESSING_ILLUMINA(reads_ch.illumina)

    emit:
    preprocessed_reads_ch = PREPROCESSING_ILLUMINA.out.preprocessed_reads_ch.mix(not_preprocessed_ch)
    collated_trimming_stats_ch = PREPROCESSING_ILLUMINA.out.collated_trimming_stats_ch
    collated_host_reads_stats_ch = PREPROCESSING_ILLUMINA.out.collated_host_reads_stats_ch
}