#!/usr/bin/env nextflow

//
// MODULES
//
include { BOWTIE2; BOWTIE2_INDEX } from './modules/bowtie2'
include { BWA; BWA_INDEX } from './modules/bwa'
include { CONVERT_TO_BAM; SAMTOOLS_SORT; INDEX_REF; INDEX_BAM as INDEX_SORTED_BAM; INDEX_BAM as INDEX_DEDUP_BAM; SAMTOOLS_STATS } from './modules/samtools'
include { BCFTOOLS_CALL; BCFTOOLS_MPILEUP; BCFTOOLS_FILTERING; BCFTOOLS_EXTRACT; PUBLISH_VCF } from './modules/bcftools'
include { PICARD_MARKDUP } from './modules/picard'
include { CURATE_CONSENSUS } from './modules/curate'
include { BAM_COVERAGE } from './modules/deeptools'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow STRAIN_MAPPER {

    take:
    ch_reads        // tuple( meta, read_1, read_2 )
    reference       // file: given reference

    main:

    //
    //BOWTIE2 WORKFLOW
    //

    if (params.mapper == "bowtie2") {
        //BOWTIE2 INDEX
        bt2_index_files = file("${reference}.bt2")
        if (bt2_index_files.isFile()) {
            Channel.fromPath(bt2_index_files)
            | collect
            | set { ch_bt2_index }

        } else {
            BOWTIE2_INDEX( reference )
            | set { ch_bt2_index }
        }

        //
        // MAPPING: Bowtie2
        //
        BOWTIE2 ( ch_reads, ch_bt2_index )
        | set { ch_mapped }

    } else if (params.mapper == "bwa") {
        //
        //BWA WORKFLOW
        //

        // BWA INDEX
        bwa_index_files = file("${reference}.amb")
        if (bwa_index_files.isFile()) {
            index_files = Channel.fromPath("${reference}{.amb,.ann,.bwt,.pac,.sa}")

            index_files
            | collect
            | map { collected_indexes -> [reference, collected_indexes]}
            | set { ch_bwa_index }

        } else {
            BWA_INDEX(reference)
            | set { ch_bwa_index }
        }

        //
        // MAPPING: Bwa
        //
        BWA( ch_reads, ch_bwa_index )
        | set { ch_mapped }

    } else {
        error "supplied mapper: ${params.mapper} is not currently supported"
    }


    // INDEX REF FASTA FOR DOWNSTREAM PROCESSES
    faidx_file = file("${reference}.fai")
    if (faidx_file.isFile()) {
        Channel.of( [reference, faidx_file] )
        | set { ch_ref_index }

    } else {
        INDEX_REF(reference)
        | set { ch_ref_index }
    }

    //
    // POST-PROCESSING
    //

    CONVERT_TO_BAM( ch_mapped )
    | SAMTOOLS_SORT
    | INDEX_SORTED_BAM
    | set { ch_sorted_reads }

    if (params.skip_read_deduplication){
        ch_sorted_reads
        | set { bam_index }

    } else {
        PICARD_MARKDUP(ch_sorted_reads)
        | INDEX_DEDUP_BAM
        | set { bam_index }
    }

    if (params.bigwig){
        BAM_COVERAGE(bam_index)
        BAM_COVERAGE.out.finished_ch
        | set { coverage_finished }
    } else {
        coverage_finished = Channel.value("BAM_COVERAGE not run")
    }

    if (params.samtools_stats){
        SAMTOOLS_STATS(bam_index)
        SAMTOOLS_STATS.out.finished_ch
        | set { stats_finished }
    } else {
        stats_finished = Channel.value("SAMTOOLS_STATS not run")
    }


    bam_index
    | combine(ch_ref_index)
    | BCFTOOLS_MPILEUP
    | BCFTOOLS_CALL
    | set { ch_vcf_allpos }

    if (!params.skip_filtering) {
        BCFTOOLS_FILTERING(ch_vcf_allpos )
        | set { ch_vcf_final }
        
        BCFTOOLS_EXTRACT(ch_vcf_final, Channel.value('FILTER="Het"'), Channel.value("heterozygous_sites"))
    } else{
        ch_vcf_allpos.set { ch_vcf_final }
    }

    PUBLISH_VCF( ch_vcf_final )
    
    ch_vcf_final
    | combine(ch_ref_index)
    | set { ch_vcf_and_ref }

    CURATE_CONSENSUS( ch_vcf_and_ref )
    CURATE_CONSENSUS.out.finished_ch
    | set { consensus_finished }

    if (!params.skip_cleanup) {
        ch_mapped.join(CONVERT_TO_BAM.out.mapped_reads_bam) // join all contents of "channel" together towards deletion
        | join(SAMTOOLS_SORT.out.sorted_reads)
        | join(INDEX_SORTED_BAM.out.indexed_bam)
        | join(BCFTOOLS_MPILEUP.out.mpileup_file)
        | join(BCFTOOLS_CALL.out.vcf_allpos)
        | join(coverage_finished)  // use these *_finished dummy value channels as a way to ensure waiting on completion of all branches of the workflow
        | join(stats_finished)  
        | join(consensus_finished)
        | flatten
        | filter(Path)
        | map { it.delete() }

        if (!params.skip_read_deduplication) {
            PICARD_MARKDUP.out.dedup_reads
            | join(coverage_finished)
            | join(stats_finished)  
            | join(consensus_finished)
            | flatten
            | filter(Path)
            | map { it.delete() }
        }
    }

    emit: 
    CURATE_CONSENSUS.out.curated_consensus
}

/*
========================================================================================
    THE END
========================================================================================
*/
