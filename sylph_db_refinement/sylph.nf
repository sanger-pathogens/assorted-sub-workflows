#!/usr/bin/env nextflow

/*
========================================================================================
    IMPORT MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULES
//
include { SYLPH_SKETCH_DB;
          SYLPH_PROFILE;
          SYLPH_SUMMARIZE } from './modules/sylph.nf'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow SYLPH_DB_REFINEMENT {
    take:
    reads_ch // meta, read_1, read_2

    main:
    def sylph_db_ch

    if (params.assemblies) {
        def assemblies_ch = channel.fromPath(params.assemblies).first()
        SYLPH_SKETCH_DB(assemblies_ch)
        // If assemblies are provided, build a fresh Sylph DB and use that.
        sylph_db_ch = SYLPH_SKETCH_DB.out.db_sketch
    } else if (params.sylph_db_custom) {
        // Otherwise prefer a user-provided custom .syldb.
        sylph_db_ch = channel.fromPath(params.sylph_db_custom).first()
    } else {
        // Fall back to the default Sylph DB from config.
        sylph_db_ch = channel.fromPath(params.sylph_db).first()
    }

    SYLPH_PROFILE(reads_ch, sylph_db_ch)
    // Summarize across all sample reports (thresholds).
    | map { meta, report -> report }
    | collect()
    | SYLPH_SUMMARIZE

    emit:
    references = SYLPH_SUMMARIZE.out.references
    sylph_summary = SYLPH_SUMMARIZE.out.sylph_summary
    sylph_db = sylph_db_ch
}

/*
========================================================================================
    THE END
========================================================================================
*/
