#!/usr/bin/env nextflow

/*
========================================================================================
    IMPORT MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULES
//
include { SYLPH_SKETCH;
          SYLPH_PROFILE;
          SYLPH_SUMMARIZE;
          SYLPHTAX_TAXPROF } from '../taxo_profile/modules/sylph.nf'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow SYLPH_REFSET {
    take:
    reads_ch // meta, read_1, read_2

    main:
    def sylph_db_ch = channel.fromPath(params.sylph_db).first()
    def sylph_tax_metadata_ch = channel.fromPath(params.sylph_tax_metadata).first()

    SYLPH_SKETCH(reads_ch)
    SYLPH_PROFILE(SYLPH_SKETCH.out.sketch)

    // Taxonomic annotations from Sylph profile output.
    SYLPHTAX_TAXPROF(SYLPH_PROFILE.out.sylph_report, sylph_tax_metadata_ch)
    def taxprof_done = SYLPHTAX_TAXPROF.out.sylphtax_mpa_report
        .map { meta, report -> true }
        .collect()
        .first()

    // Summarize across all sample reports (thresholds).
    SYLPH_PROFILE.out.sylph_report
        .map { meta, report -> report }
        .collect()
        .combine(taxprof_done)
        .map { reports, done -> reports }
        | SYLPH_SUMMARIZE

    emit:
    references = SYLPH_SUMMARIZE.out.references
    sylph_summary = SYLPH_SUMMARIZE.out.sylph_summary
    sylph_db = sylph_db_ch
    sylphtax_mpa_report = SYLPHTAX_TAXPROF.out.sylphtax_mpa_report
}

/*
========================================================================================
    THE END
========================================================================================
*/
