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
            SYLPHTAX_TAXPROF } from '../taxo_profile/modules/sylph.nf'
include   { COMBINE_SYLPH_REPORTS;
            NORMALIZE_QUERY_REPORT_FOR_SYLPHTAX;
            SYLPH_FILTER;
            GROUP_SYLPH_REFS_BY_TAXON;
            COMBINE_REFS_ACROSS_SAMPLES
            EXPAND_REFS } from './modules/helper_processes.nf'

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

    SYLPH_SKETCH(reads_ch)
    SYLPH_QUERY(SYLPH_SKETCH.out.sketch)
    sylph_report_ch = SYLPH_QUERY.out.sylph_report

    // Combine sylph reports
    sylph_report_ch
    | map { meta, report -> report }
    | collectFile( name: "combined_sylph_report.tsv", keepHeader: true, storeDir: "${params.outdir}/sylph" )
    | map { report -> [ [ID: "combined"], report ] }
    | set { combined_sylph_report }

    // Filter based on ANI and coverage threshold
    combined_sylph_report
    | combine(sylph_tax_metadata_ch)
    | map { meta, report, taxonomy_data -> [ meta, report, taxonomy_data ] }
    | SYLPH_FILTER

    // Normalize query output to make it compatible for input to sylph-tax.
    NORMALIZE_QUERY_REPORT_FOR_SYLPHTAX(SYLPH_FILTER.out.report)
    sylphtax_input_ch = NORMALIZE_QUERY_REPORT_FOR_SYLPHTAX.out.sylph_report

    // Get taxonomic profile in metaphlan (mpa) report format
    sylphtax_input_ch
    | combine(sylph_tax_metadata_ch)
    | SYLPHTAX_TAXPROF

    // Get all references for taxonomic groups output from sylph-tax report
    if (params.expand_refs) {
        genome_id_to_file = channel.fromPath(params.genome_id_to_file).first()

        SYLPHTAX_TAXPROF.out.sylphtax_mpa_report
        | combine(sylph_tax_metadata_ch)
        | combine(genome_id_to_file)
        | EXPAND_REFS

        EXPAND_REFS.out.references
        | transpose
        | map { meta, ref ->
            def new_meta = ["ID": ref.baseName]  // Construct taxonomic group from filename
            [ new_meta, ref ]
        }
        | set { references }

    } else {
        // Join sylph report (incl. reference filepaths) with taxonomy
        sylphtax_input_ch
        | join(SYLPHTAX_TAXPROF.out.sylphtax_mpa_report)
        | GROUP_SYLPH_REFS_BY_TAXON

        // Group reports by taxonomic group and output taxon-specific references
        GROUP_SYLPH_REFS_BY_TAXON.out.taxon_group_refs
        | transpose
        | map { meta, ref ->
            def new_meta = ["ID": ref.baseName]  // Construct taxonomic group from filename
            [new_meta, ref]
        }
        | set { references }

    }

    emit:
    references
    combined_sylph_report = SYLPH_FILTER.out.report
    //TODO Do we need to output the following?
    // sylph_db = sylph_db_ch
    // sylphtax_mpa_report = SYLPHTAX_TAXPROF.out.sylphtax_mpa_report
}

/*
========================================================================================
    THE END
========================================================================================
*/
