#!/usr/bin/env nextflow

/*
========================================================================================
    IMPORT MODULES/SUBWORKFLOWS
========================================================================================
*/


//
// MODULES
//
include { SYLPH_SKETCH 
          SYLPH_PROFILE
          SYLPHTAX_TAXPROF } from './modules/sylph.nf'
include { CALC_DIVERSITY    } from './modules/diversity.nf'

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
    def sylph_tax_metadata_ch = channel.fromPath(params.sylph_tax_metadata).first()

    if (params.sketch_size) {
        log.warn "Parameter '--sketch_size' is deprecated. Use '--sylph_k' instead."
        params.sylph_k = params.sketch_size
    }
    
    if (params.sylph_profile){
        SYLPH_SKETCH(read_ch)
        | SYLPH_PROFILE
        | set { sylph_report }

        sylph_report
        | combine(sylph_tax_metadata_ch)
        | SYLPHTAX_TAXPROF
        | set { sylphtax_mpa_report }
    } else {
        sylph_report = channel.empty()
        sylphtax_mpa_report = channel.empty()
    }
    if (params.bracken_profile){
        KRAKEN2BRACKEN(read_ch)
        ch_kraken2_style_bracken_reports = KRAKEN2BRACKEN.out.ch_kraken2_style_bracken_reports
        ch_mpa_abundance_reports = KRAKEN2BRACKEN.out.ch_mpa_abundance_reports
    } else {
        ch_kraken2_style_bracken_reports = channel.empty()
        ch_mpa_abundance_reports = channel.empty()
    }
    if (params.calculate_alpha_diverity) {
        sylph_report
        | map { meta, file -> file.toString() }   
        | collectFile(name: 'sylph_profiles.txt', newLine: true)
        | set{ diversity_in }
        CALC_DIVERSITY(diversity_in)      
    }

    emit:
    sylph_report
    sylphtax_mpa_report
    ch_kraken2_style_bracken_reports
    ch_mpa_abundance_reports
}

/*
========================================================================================
    THE END
========================================================================================
*/
