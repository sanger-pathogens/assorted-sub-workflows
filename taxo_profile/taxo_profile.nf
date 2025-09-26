#!/usr/bin/env nextflow

/*
========================================================================================
    IMPORT MODULES/SUBWORKFLOWS
========================================================================================
*/


//
// MODULES
//
include { SYLPH_PROFILE; SYLPH_SKETCH } from './modules/sylph.nf'

//
// SUBWORKFLOWS
//
include { KRAKEN2BRACKEN } from '../kraken2bracken/subworkflows/kraken2bracken.nf'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow TAXO_PROFILE {
    take:
    read_ch // meta, read_1, read_2

    main:
    if (params.sylph_profile){
        SYLPH_SKETCH(read_ch)
        | SYLPH_PROFILE
        | set { sylph_report }
    } else {
        sylph_report = channel.empty()
    }
    if (params.bracken_profile){
        KRAKEN2BRACKEN(read_ch)
        ch_kraken2_style_bracken_reports = KRAKEN2BRACKEN.out.ch_kraken2_style_bracken_reports
    } else {
        ch_kraken2_style_bracken_reports = channel.empty()
    }

    emit:
    sylph_report
    ch_kraken2_style_bracken_reports
}

/*
========================================================================================
    THE END
========================================================================================
*/
