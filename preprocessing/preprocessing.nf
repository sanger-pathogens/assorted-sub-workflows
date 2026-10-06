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
    // but it is a safety measure to avoid processing the wring type of data. 
    // We could drop this in the future when settled on how we handle mutliple patforms.
    reads_ch.branch{ meta, fwd, rev ->
        illumina: meta.Platform == "ILLUMINA"

        ONT: meta.Platform == "ONT"

        other: true
    }
    | set { reads_ch }

    PREPROCESSING_ILLUMINA(reads_ch.illumina) // in the future, can insert here other preprocessing workflows for ONT or other platforms

    emit:
    // in the future, we may want to combine these output channels with those derived from an PREPROCESSING_ONT workflow here,
    // but for now, just emmit the PREPROCESSING_ILLUMINA outputs directly
    preprocessed_reads_ch = PREPROCESSING_ILLUMINA.out.preprocessed_reads_ch
    collated_trimming_stats_ch = PREPROCESSING_ILLUMINA.out.collated_trimming_stats_ch
    collated_host_reads_stats_ch = PREPROCESSING_ILLUMINA.out.collated_host_reads_stats_ch
}