#!/usr/bin/env nextflow

/*
========================================================================================
    IMPORT MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULES
//
include { KRAKEN2; KRAKEN2_GET_CLASSIFIED; COMPRESS_READS } from '../modules/kraken2bracken/kraken2'
include { BRACKEN } from '../modules/kraken2bracken/bracken'
include { KREPORT2MPA; GENERATE_ABUNDANCE_SUMMARY } from '../modules/kraken2bracken/krakentools'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow KRAKEN2BRACKEN{

    take:
    ch_reads // meta, paired_reads

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
    // Assume a pre-built bracken database file has been generated from the given kraken2 database and moved into this database directory
    required_kmer_distrib = file("${kraken2_db_dir}/database${params.read_len}mers.kmer_distrib")
    if (!required_kmer_distrib.exists()) {
        log.error("The required bracken kmer distribution database file cannot be found in the kraken database directory ${kraken2_db_dir}")
        exit 1
    } else {
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
        .dump(tag: 'mpa_abundance_reports')
        .set { ch_mpa_abundance_reports }
    GENERATE_ABUNDANCE_SUMMARY(
        ch_mpa_abundance_reports
    )

    emit:
    BRACKEN.out.kraken_style_bracken_report
}

/*
========================================================================================
    THE END
========================================================================================
*/
