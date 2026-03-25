#!/usr/bin/env nextflow

/*
========================================================================================
    IMPORT MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULES
//
include   { SYLPH_SKETCH;
            SYLPH_QUERY;
            SYLPH_PROFILE;
            SYLPH_SUMMARIZE;
            SYLPHTAX_TAXPROF } from '../taxo_profile/modules/sylph.nf'
include   { COMBINE_SYLPH_REPORTS;
            NORMALIZE_QUERY_REPORT_FOR_SYLPHTAX;
            GROUP_SYLPH_REFS_BY_TAXON;
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

    sylph_method = (params.sylph_method ?: 'query').toString().toLowerCase()

    SYLPH_SKETCH(reads_ch)
    if (sylph_method == 'query') {
        SYLPH_QUERY(SYLPH_SKETCH.out.sketch)
        sylph_report_ch = SYLPH_QUERY.out.sylph_report
    } else if (sylph_method == 'profile') {
        SYLPH_PROFILE(SYLPH_SKETCH.out.sketch)
        sylph_report_ch = SYLPH_PROFILE.out.sylph_report
    } else {
        error "Unsupported --sylph_method '${params.sylph_method}'. Use 'query' or 'profile'."
    }


    sylph_report_ch
    | map { meta, report -> report }
    | collect()
    | map { reports -> [[ID: 'all_samples'], reports] }
    | COMBINE_SYLPH_REPORTS

    if (sylph_method == 'query') {
        NORMALIZE_QUERY_REPORT_FOR_SYLPHTAX(COMBINE_SYLPH_REPORTS.out.sylph_report)
        sylphtax_input_ch = NORMALIZE_QUERY_REPORT_FOR_SYLPHTAX.out.sylph_report
    } else {
        sylphtax_input_ch = COMBINE_SYLPH_REPORTS.out.sylph_report
    }

    // Combine sylph reports
    SYLPH_QUERY.out.sylph_report
    | map { meta, report -> report }
    | collectFile( name: "combined_sylph_report.tsv", keepHeader: true )
    // | map { reports -> [ "combined_sylph_report", (reports instanceof List) ? reports : [reports] ] }
    | map { report -> [ [ID: "combined_sylph_report"], report ] }
    | set { combined_sylph_report }
    
    // COMBINE_SYLPH_REPORTS(
    //     combined_sylph_report,
    //     "${params.outdir}/sylph/"
    // )
    
    // Filter based on ANI and coverage threshold
    // COMBINE_SYLPH_REPORTS.out.group_report
    // | map { group, report ->
    //     def meta = [:]
    //     meta.ID = group
    //     [ meta, report ]
    // }
    combined_sylph_report
    | SYLPH_SUMMARIZE

    // Get taxonomic profile in metaphlan (mpa) report format
    // COMBINE_SYLPH_REPORTS.out.group_report
    SYLPH_SUMMARIZE.out.report
    | NORMALIZE_SYLPH_QUERY_REPORT
    | combine(sylph_tax_metadata_ch)
    | SYLPHTAX_TAXPROF

    // Join sylph report (incl. reference filepaths) with taxonomy
    SYLPH_SUMMARIZE.out.report
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
    | groupTuple
    | set { grouped_sample_reports }

    // Extract reference lists from reports?

    // COMBINE_REFS_ACROSS_SAMPLES(
    //     grouped_sample_reports,
    //     null // publish_dir (don't publish)
    // )
    // | map { taxon_group, report ->
    //     def meta = [:]
    //     meta.ID = taxon_group
    //     [meta, report]
    // }
    // | SYLPH_SUMMARIZE

    // SYLPH_SUMMARIZE.out.sylph_summary
    // | map { meta, summary -> summary }
    // | collect
    // | map { summaries ->
    //     [ "summaries", summaries ]
    // }
    // | set { sylph_summaries }

    // COMBINE_SYLPH_SUMMARIES(
    //     sylph_summaries,
    //     "${params.outdir}/sylph/summaries"  // publishDir
    // )
    // COMBINE_SYLPH_REFERENCES(
    //     ,
    // )

    emit:
    references = SYLPH_SUMMARIZE.out.references
    sylph_summary = SYLPH_SUMMARIZE.out.sylph_summary
    combined_sylph_report = COMBINE_SYLPH_REPORTS.out.sylph_report
    //TODO Do we need to output the following?
    // sylph_db = sylph_db_ch
    // sylphtax_mpa_report = SYLPHTAX_TAXPROF.out.sylphtax_mpa_report
}

/*
========================================================================================
    THE END
========================================================================================
*/
