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
include { KRAKEN2BRACKEN } from '../kraken2bracken/subworkflows/kraken2bracken.nf'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow QC {
    take:
    ch_reads // meta, read_1, read_2

    main:
    reads_ch
    | FASTQC
    | PASS_OR_FAIL_FASTQC
    | set { fastqc_results }

    reads_ch
    | KRAKEN2BRACKEN
    | PASS_OR_FAIL_K2B
    | set { k2b_results }

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
