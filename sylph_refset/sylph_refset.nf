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

workflow SYLPH_REF_SELECTION {
    take:
    reads_ch // meta, read_1, read_2

    main:
    sylph_db_ch = channel.fromPath(params.sylph_db).first()
    sylph_tax_metadata_ch = channel.fromPath(params.sylph_tax_metadata).first()

    SYLPH_SKETCH(reads_ch)
    | SYLPH_QUERY
    
    SYLPH_QUERY.out.sylph_report
        .combine(sylph_tax_metadata_ch)
        | SYLPHTAX_TAXPROF

    SYLPHTAX_TAXPROF.out.sylphtax_mpa_report
        .map { true }
        .collect()
        .first()
        | set { taxprof_done }

    SYLPH_QUERY.out.sylph_report
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
