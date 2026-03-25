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
include { GROUP_SYLPH_REFS_BY_TAXON;
          COMBINE_REFS_ACROSS_SAMPLES } from './modules/helper_processes.nf'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow SYLPH_REF_SELECTION {
    take:
    reads_ch // meta, read_1, read_2

    main:
    // Workflow only works with GTDB-like sylph database
    sylph_db_ch = channel.fromPath(params.sylph_db).first()
    sylph_tax_metadata_ch = channel.fromPath(params.sylph_tax_metadata).first()

    // Profile sample using sylph
    SYLPH_SKETCH(reads_ch)
    | SYLPH_QUERY

    // Get taxonomic profile in metaphlan (mpa) report format
    SYLPH_QUERY.out.sylph_report
    | combine(sylph_tax_metadata_ch)
    | SYLPHTAX_TAXPROF

    // Join sylph report (incl. reference filepaths) with taxonomy
    SYLPH_QUERY.out.sylph_report
    | join(SYLPHTAX_TAXPROF.out.sylphtax_mpa_report)
    | GROUP_SYLPH_REFS_BY_TAXON

    // Group reports by taxonomic group, then combine them across samples
    // Output can then be passed to SYLPH_SUMMARIZE for filtering
    GROUP_SYLPH_REFS_BY_TAXON.out.taxon_group_ref_reports
    | transpose
    | map { meta, report ->
        taxon_group = report.baseName.replace("${meta.ID}_", "")  // Construct taxonomic group from filename
        [taxon_group, report]
    }
    | groupTuple()
    | COMBINE_REFS_ACROSS_SAMPLES

    COMBINE_REFS_ACROSS_SAMPLES.out.taxon_group_ref_report
    | map { taxon_group, report ->
        def meta = [:]
        meta.ID = taxon_group
        [meta, report]
    }
    | SYLPH_SUMMARIZE

    emit:
    references = SYLPH_SUMMARIZE.out.references
    sylph_summary = SYLPH_SUMMARIZE.out.sylph_summary
    //TODO Do we need to output the following?
    // sylph_db = sylph_db_ch
    // sylphtax_mpa_report = SYLPHTAX_TAXPROF.out.sylphtax_mpa_report
}

/*
========================================================================================
    THE END
========================================================================================
*/
