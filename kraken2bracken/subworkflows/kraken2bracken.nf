#!/usr/bin/env nextflow

/*
========================================================================================
    IMPORT MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULES
//
include { KRAKEN2; KRAKEN2_GET_CLASSIFIED; COMPRESS_READS } from '../modules/kraken2'
include { BRACKEN } from '../modules/bracken'
include { KREPORT2MPA; GENERATE_ABUNDANCE_SUMMARY } from '../modules/krakentools'

//
// SUBWORKFLOWS
//
include { BRACKEN_BUILD } from './bracken_build.nf'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow KRAKEN2BRACKEN{

    take:
    ch_reads // meta, read_1, read_2

    main:

    Channel.fromPath(params.kraken2_db)
        .set { ch_kraken2_db }

    //
    // CLASSIFICATION
    //
    ch_reads
        .combine(ch_kraken2_db)
        .dump(tag: 'reads_and_kraken2_db')
        .set { ch_reads_and_kraken2_db }
    
    if (params.get_classified_reads) {
        KRAKEN2_GET_CLASSIFIED(ch_reads_and_kraken2_db)

        KRAKEN2_GET_CLASSIFIED.out.kraken2_sample_report.dump(tag: 'kraken2_sample_report').set { ch_kraken2_sample_report }
        KRAKEN2_GET_CLASSIFIED.out.classified_reads.dump(tag: 'kraken2_classified_reads').set { ch_kraken2_classified_reads }
        KRAKEN2_GET_CLASSIFIED.out.unclassified_reads.dump(tag: 'kraken2_unclassified_reads').set { ch_kraken2_unclassified_reads }
        ch_kraken2_classified_reads.join(ch_kraken2_unclassified_reads)
            .map { meta, classified, unclassified -> [meta, classified + unclassified] }
            .set { ch_reads_to_compress }

        COMPRESS_READS(
            ch_reads_to_compress
        )
    } else {
        KRAKEN2(ch_reads_and_kraken2_db)
        KRAKEN2.out.kraken2_sample_report.dump(tag: 'kraken2_sample_report').set { ch_kraken2_sample_report }
    }

    //
    // ABUNDANCE ESTIMATION
    //

    kraken2_db_dir = file(params.kraken2_db, checkIfExists: true)

    // Define the required k-mer distribution file path
    required_kmer_distrib = file("${kraken2_db_dir}/database${params.read_len}mers.kmer_distrib")

    // Check if building is enabled and the required file doesn't exist
    if (params.enable_building && !required_kmer_distrib.exists()) {
        BRACKEN_BUILD(ch_kraken2_db)
        BRACKEN_BUILD.out.ch_kmer_distrib
            .dump(tag: 'kmer_distrib')
            .set { ch_kmer_distrib }
    }
    // If building is not enabled and the file doesn't exist, log an error
    else if (!params.enable_building && !required_kmer_distrib.exists()) {
        error("Required k-mer distribution file does not exist, and building is not enabled. Please enable building with --enable_building.")
    }
    // If the required file exists, load it into a channel
    else {
        Channel.fromPath(required_kmer_distrib)
            .dump(tag: 'kmer_distrib')
            .set { ch_kmer_distrib }
    }

    ch_kraken2_sample_report
        .combine(ch_kmer_distrib)
        .dump(tag: 'kraken2_report_and_kmer_distrib')
        .set { ch_kraken2_report_and_kmer_distrib }
    BRACKEN(
        ch_kraken2_report_and_kmer_distrib
    )

    KREPORT2MPA(BRACKEN.out.kraken_style_bracken_report)

    //
    // SUMMARISE ABUNDANCE
    //
    KREPORT2MPA.out.mpa_abundance_report
        .map { meta, report -> report }
        .collect()
        .unique()
        .dump(tag: 'mpa_abundance_reports')
        .set { ch_mpa_abundance_reports }
    GENERATE_ABUNDANCE_SUMMARY(
        ch_mpa_abundance_reports
    )

    emit:
    ch_kraken2_style_bracken_reports = BRACKEN.out.kraken_style_bracken_report
    kraken2_report_for_multiqc = ch_kraken2_report_and_kmer_distrib
}

/*
========================================================================================
    THE END
========================================================================================
*/
