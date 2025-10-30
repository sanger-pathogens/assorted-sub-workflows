#!/usr/bin/env nextflow

include { TRIMMOMATIC        } from "./modules/trimmomatic.nf"
include { FASTQ2FASTA        } from "./modules/fastq2fasta.nf"
include { TRF                } from "./modules/trf.nf"
include { RMREPEATFROMFASTQ  } from "./modules/rmRepeatFromFq.nf"
include { BMTAGGER           } from "./modules/bmtagger.nf"
include { FILTER_HOST_READS; 
          GET_HOST_READS     } from './modules/filter_reads.nf'
include { GENERATE_STATS     } from './modules/generate_stats.nf'
include { COLLATE_STATS      } from './modules/collate_stats.nf'
include { COMPRESS_READS; 
          RENAME_READS       } from "./modules/helper_processes.nf"
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
    if (params.run_trimmomatic){
        TRIMMOMATIC(reads_ch)
        trimm_Out_ch = TRIMMOMATIC.out.results // tuple (meta, read_trim_1, read_trim_2)
    }

    if (params.run_trf){
        if (params.run_trimmomatic){
            fastq2fasta_in_ch = trimm_Out_ch
        } else {
            fastq2fasta_in_ch = reads_ch
        }
        
        fastq2fasta_in_ch.map{
            meta, fq_1, fq_2 -> 
            tuple (meta.id, meta, [fq_1, fq_2])
        }
        .set {fqs_ch}
        
        FASTQ2FASTA(fastq2fasta_in_ch)
        
        TRF(FASTQ2FASTA.out) // tuple (meta, fasta_1, fasta_2)
        TRF.out.results | map {
            meta, trf_out_1, trf_out_2 -> 
            tuple (meta.id,[trf_out_1, trf_out_2])
        }
        | set {trf_Out_ch}  // tuple (meta.id, [trfs_out])

        fqs_ch // tuple (meta.id, meta, [fqs])
        | join(trf_Out_ch) // tuple (meta.id, meta, [fqs], [trfs_out])
        | map {id, meta, fqs, trfs -> tuple(meta, fqs[0], fqs[1], trfs[0], trfs[1])}
        | set {rmRepeat_In_ch}

        RMREPEATFROMFASTQ(rmRepeat_In_ch)
        RMREPEATFROMFASTQ.out.fastqs
            | set {rmRepeat_Out_ch} // tuple (meta, trf_fq_1, trf_fq_2)
    }


    // run bmtagger
    if (params.run_hrr){

        // if only scrubber is on
        if ((!params.run_trimmomatic) && (!params.run_trf)){
            bmtagger_In_ch = reads_ch
        }

        // if trimmomatic and no trf
        if ((params.run_trimmomatic) && (!params.run_trf)){
            bmtagger_In_ch = trimmomatic_Out_ch
        }

        // if trf is true
        if (params.run_trf){
            bmtagger_In_ch = trf_Out_ch
        }

            // run hrr
        BMTAGGER(bmtagger_In_ch)
        | set {bmtagger_Out_ch}

        FILTER_HOST_READS(bmtagger_Out_ch.data_ch, bmtagger_Out_ch.bmtagger_list_ch)
        | set { filtered_reads_Out_ch }

        GET_HOST_READS(bmtagger_Out_ch.data_ch, bmtagger_Out_ch.bmtagger_list_ch)
        | set { host_reads_Out_ch }

        filtered_reads_Out_ch.data_ch
        | join(filtered_reads_Out_ch.cleaned_ch)
        | join(host_reads_Out_ch.host_ch)
        | join(reads_ch)
        | GENERATE_STATS

        COLLATE_STATS(GENERATE_STATS.out.stats_ch.collect())
    }

    // setup output channel
    // if trimmomatics on, trf off and scrubber off
    if ((params.run_trimmomatic) && (!params.run_trf) && (!params.run_hrr)){
        cleaned_reads_ch = trimmomatic_Out_ch // tuple (meta, reads_trim_1, reads_trim_2)
    }
    // if trf on, scrubber off
    if ((params.run_trf) && (!params.run_hrr)){
        cleaned_reads_ch = trf_Out_ch // tuple (meta, trf_fq_1, trf_fq_2)
    }
    // if scrubber on
    if (params.run_hrr){
        cleaned_reads_ch = filtered_reads_Out_ch.cleaned_ch // tuple (meta, reads_clean_1, reads_clean_2)
    }

    // publish compressed clean reads
    if (params.publish_clean_reads){
        if (params.compress_clean_reads){
            COMPRESS_READS(cleaned_reads_ch)
        } else {
            RENAME_READS(cleaned_reads_ch)
        }
    }
    
    FASTQC(cleaned_reads_ch) 

    fastqc_pass_criteria = file(params.fastqc_pass_criteria, checkIfExists: true)
    fastqc_no_fail_criteria = file(params.fastqc_no_fail_criteria, checkIfExists: true)

    PASS_OR_FAIL_FASTQC(FASTQC.out.zip, fastqc_pass_criteria, fastqc_no_fail_criteria)
    | set { fastqc_results }

    TAXO_PROFILE(cleaned_reads_ch)

    FASTQC.out.zip.collect{it[1,2]}
        | mix(TAXO_PROFILE.out.ch_kraken2_style_bracken_reports.collect{it[1]})
        | collect
        | set { multiqc_input }

    pass_fail_channel = fastqc_results.map { sample_pass_fail -> [sample_pass_fail[0].ID, sample_pass_fail[1]] } // Sample ID, Pass/Fail status e.g [ERR14241855, pass] 

    if (params.bracken_profile) {

        TAXO_PROFILE.out.ch_kraken2_style_bracken_reports
        | PASS_OR_FAIL_K2B
        | set { k2b_results }

        pass_fail_channel = pass_fail_channel.join(k2b_results.map { sample_pass_fail -> [sample_pass_fail[0].ID, sample_pass_fail[1]] })
    
    }

    if (params.sylph_profile) {
        
        TAXO_PROFILE.out.sylphtax_mpa_report
            | PASS_OR_FAIL_SYLPH
            | set { sylph_results }

        pass_fail_channel = pass_fail_channel.join(sylph_results.map { sample_pass_fail -> [sample_pass_fail[0].ID, sample_pass_fail[1]] })

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
