#!/usr/bin/env nextflow
include { PREPROCESSING  } from '../preprocessing/preprocessing.nf'
include { FASTQC              } from './modules/fastqc.nf'
include { PASS_OR_FAIL_FASTQC
          PASS_OR_FAIL_K2B
          PASS_OR_FAIL_SYLPH  } from './modules/pass_or_fail.nf'
include { REPORT              } from './modules/reporting.nf'
include { TAXO_PROFILE        } from '../taxo_profile/taxo_profile.nf'

workflow QC {
    take:
    reads_ch // meta, read_1, read_2

    main:

    if (!params.skip_preprocessing) {
        PREPROCESSING(reads_ch)
        | FASTQC 
    }
    else {
        FASTQC(reads_ch)
    }
    
    

    fastqc_pass_criteria = file(params.fastqc_pass_criteria, checkIfExists: true)
    fastqc_no_fail_criteria = file(params.fastqc_no_fail_criteria, checkIfExists: true)

    PASS_OR_FAIL_FASTQC(FASTQC.out.zip, fastqc_pass_criteria, fastqc_no_fail_criteria)
    | set { fastqc_results }

    TAXO_PROFILE(reads_ch)

    FASTQC.out.zip.collect{it[1,2]}
        | mix(TAXO_PROFILE.out.ch_kraken2_style_bracken_reports.collect{it[1]})
        | collect
        | set { multiqc_input }

    pass_fail_fastqc_ch = fastqc_results.map { sample_pass_fail -> [sample_pass_fail[0].ID, sample_pass_fail[1]] } // Sample ID, Pass/Fail status e.g [ERR14241855, pass] 

    if (params.bracken_profile) {

        TAXO_PROFILE.out.ch_kraken2_style_bracken_reports
        | PASS_OR_FAIL_K2B
        | set { k2b_results }

        pass_fail_k2b_ch = k2b_results.map { sample_pass_fail -> [sample_pass_fail[0].ID, sample_pass_fail[1]] }
    }

    if (params.sylph_profile) {
        
        TAXO_PROFILE.out.sylphtax_mpa_report
            | PASS_OR_FAIL_SYLPH
            | set { sylph_results }

        pass_fail_sylph_ch = sylph_results.map { sample_pass_fail -> [sample_pass_fail[0].ID, sample_pass_fail[1]] }

    }

    if (params.bracken_profile && params.sylph_profile) {
        pass_fail_channel = pass_fail_fastqc_ch.join(pass_fail_k2b_ch)
                            .map { flat -> flat.flatten() }
                            .join(pass_fail_sylph_ch)
                            .map { flat -> flat.flatten() }
    
    } 
    else if (params.bracken_profile) {
        pass_fail_channel = pass_fail_fastqc_ch.join(pass_fail_k2b_ch)
                            .map { flat -> flat.flatten() }
    } 
    else if (params.sylph_profile) {
        pass_fail_channel = pass_fail_fastqc_ch.join(pass_fail_sylph_ch)
                            .map { flat -> flat.flatten() }
    } 
    else {
        pass_fail_channel = pass_fail_fastqc_ch
    }

    // map dynamic number of columns into a list of sample lists for reporting
    pass_fail_channel = pass_fail_channel.collect()
                        .map { flat -> 
                        int numCols = 2 + (params.bracken_profile ? 1 : 0) + (params.sylph_profile ? 1 : 0)
                        flat.collate(numCols) }
    
    REPORT(pass_fail_channel)

    emit:
    multiqc_input
    kraken2_style_bracken_reports = TAXO_PROFILE.out.ch_kraken2_style_bracken_reports
    bracken_mpa_reports = TAXO_PROFILE.out.ch_mpa_abundance_reports
    sylphtax_mpa_report = TAXO_PROFILE.out.sylphtax_mpa_report
    qc_summary = REPORT.out.qc_summary
}
