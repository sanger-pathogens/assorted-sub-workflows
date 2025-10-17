#!/usr/bin/env nextflow

/*
========================================================================================
    IMPORT MODULES/SUBWORKFLOWS
========================================================================================
*/


//
// MODULES
//
include { FASTQC } from './modules/fastqc.nf'
include { PASS_OR_FAIL_FASTQC; PASS_OR_FAIL_K2B } from './modules/pass_or_fail.nf'

//
// SUBWORKFLOWS
//
include { TAXO_PROFILE } from '../taxo_profile/taxo_profile.nf'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow QC {
    take:
    read_ch // meta, read_1, read_2

    main:
    FASTQC(read_ch)

    fastqc_pass_criteria = file(params.fastqc_pass_criteria, checkIfExists: true)
    fastqc_no_fail_criteria = file(params.fastqc_no_fail_criteria, checkIfExists: true)

    PASS_OR_FAIL_FASTQC(FASTQC.out.zip, fastqc_pass_criteria, fastqc_no_fail_criteria)
    | set { fastqc_results }

    TAXO_PROFILE(read_ch)

    TAXO_PROFILE.out.ch_kraken2_style_bracken_reports
    | PASS_OR_FAIL_K2B
    | set { k2b_results }

    TAXO_PROFILE.out.sylphtax_mpa_report
    | PASS_OR_FAIL_SYLPH
    | set { sylph_results }

    fastqc_results
    | join(k2b_results)
    | set { results }

    emit:
    results
}

/*
========================================================================================
    THE END
========================================================================================
*/
