#!/usr/bin/env nextflow

/*
========================================================================================
    IMPORT MODULES/SUBWORKFLOWS
========================================================================================
*/


//
// MODULES
//
include { FASTQC              } from './modules/fastqc.nf'
include { PASS_OR_FAIL_FASTQC
          PASS_OR_FAIL_K2B
          PASS_OR_FAIL_SYLPH  } from './modules/pass_or_fail.nf'
include { REPORT              } from './modules/reporting.nf'

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

    FASTQC.out.zip.collect{it[1,2]}
        | mix(TAXO_PROFILE.out.ch_kraken2_style_bracken_reports.collect{it[1]})
        | collect
        | set { multiqc_input }

    pass_fail_channel = fastqc_results.map { entry -> [entry[0].ID, entry[1]] }

    if (params.bracken_profile) {

        TAXO_PROFILE.out.ch_kraken2_style_bracken_reports
        | PASS_OR_FAIL_K2B
        | set { k2b_results }

        pass_fail_channel = pass_fail_channel.join(k2b_results.map { entry -> [entry[0].ID, entry[1]] })
    
    }

    if (params.sylph_profile) {
        
        TAXO_PROFILE.out.sylphtax_mpa_report
            | PASS_OR_FAIL_SYLPH
            | set { sylph_results }

        pass_fail_channel = pass_fail_channel.join(sylph_results.map { entry -> [entry[0].ID, entry[1]] })

    }

    // map dynamic number of columns into a list of sample lists for reporting
    pass_fail_channel = pass_fail_channel.collect()
                        .map { flat -> 
                        int numCols = 2 + (params.bracken_profile ? 1 : 0) + (params.sylph_profile ? 1 : 0)
                        flat.collate(numCols) }
    
    REPORT(pass_fail_channel)

    emit:
    multiqc_input
}

/*
========================================================================================
    THE END
========================================================================================
*/
