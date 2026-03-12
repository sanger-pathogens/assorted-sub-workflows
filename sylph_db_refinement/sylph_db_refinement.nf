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
          SYLPH_QUERY;
          SYLPH_SUMMARIZE;
          SYLPHTAX_TAXPROF } from '../taxo_profile/modules/sylph.nf'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow SYLPH_DB_REFINEMENT {
    take:
    reads_ch // meta, read_1, read_2

    main:
    def sylph_db_ch = channel.fromPath(params.sylph_db).first()

    SYLPH_SKETCH(reads_ch)
    SYLPH_QUERY(SYLPH_SKETCH.out.sketch)

    // Taxonomic annotations from Sylph query output.
    SYLPHTAX_TAXPROF(SYLPH_QUERY.out.sylph_report)

    // Summarize across all sample reports (thresholds).
    SYLPH_QUERY.out.sylph_report
    | map { meta, report -> report }
    | collect()
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
