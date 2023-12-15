#!/usr/bin/env nextflow


/*
========================================================================================
    REQUIRED PARAMS IN NEXTFLOW.CONFIG
========================================================================================
*/

//   only_report_alts = true
//   VCF_filters = 'QUAL>=50 & MIN(DP)>=8 & ((ALT!="." & DP4[2]>3 & DP4[3]>3) | (ALT="." & DP4[0]>3 & DP4[1]>3))'
//   skip_filtering = false
//   keep_raw_vcf = false

//
// MODULES
//
include { BOWTIE2; BOWTIE2_INDEX } from './modules/bowtie2'
include { CONVERT_TO_BAM; SAMTOOLS_SORT; INDEX_REF } from './modules/samtools'
include { BCFTOOLS_CALL; BCFTOOLS_MPILEUP; BCFTOOLS_FILTERING; FINAL_VCF; RAW_VCF } from './modules/bcftools'
include { CURATE_CONSENSUS } from './modules/curate'

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
    // BOWTIE2 INDEX
    ref_without_extension = "${reference.parent}/${reference.baseName}"
    bt2_index_files = file("${ref_without_extension}*.bt2")
    if (bt2_index_files) {
        Channel.fromPath(bt2_index_files)
            .collect()
            .dump(tag: 'bt2_index')
            .set { ch_bt2_index }
    } else {
        BOWTIE2_INDEX(
            reference
        )
        BOWTIE2_INDEX.out.bt2_index.dump(tag: 'bt2_index').set { ch_bt2_index }
    }

    // INDEX REF FASTA
    faidx_file = file("${reference}.fai")
    if (faidx_file.isFile()) {
        Channel.of( [reference, faidx_file] ).dump(tag: 'ref_index').set { ch_ref_index }
    } else {
        INDEX_REF(
            reference
        )
        INDEX_REF.out.ref_index.dump(tag: 'ref_index').set { ch_ref_index }
    }
    
    //
    // MAPPING: Bowtie2
    //
    BOWTIE2 (
        ch_reads,
        ch_bt2_index 
    )
    BOWTIE2.out.mapped_reads.dump(tag: 'bowtie2').set { ch_mapped }

    //
    // POST-PROCESSING
    //
    CONVERT_TO_BAM(
        ch_mapped
    )
    CONVERT_TO_BAM.out.mapped_reads_bam.dump(tag: 'convert_to_bam').set { ch_mapped_reads_bam }

    SAMTOOLS_SORT(
        ch_mapped_reads_bam
    )
    SAMTOOLS_SORT.out.sorted_reads.dump(tag: 'sorted_reads').set { ch_sorted_reads }

    PICARD_MARKDUP(
        ch_sorted_reads
    )
    PICARD_MARKDUP.out.dedup_reads
        .combine(ch_ref_index)
        .dump(tag: 'sorted_reads_and_ref')
        .set { sorted_reads_and_ref }

    BCFTOOLS_MPILEUP(
        sorted_reads_and_ref
    )
    BCFTOOLS_MPILEUP.out.mpileup_file.dump(tag: 'mpileup_file').set { ch_mpileup_file }

    BCFTOOLS_CALL(
        ch_mpileup_file
    )
    BCFTOOLS_CALL.out.vcf_allpos.dump(tag: 'vcf_allpos').set { ch_vcf_allpos }

    if (params.keep_raw_vcf && !params.skip_filtering){
        RAW_VCF(
            ch_vcf_allpos
        )
    }

    if (!params.skip_filtering) {
        BCFTOOLS_FILTERING(
            ch_vcf_allpos
        )
        BCFTOOLS_FILTERING.out.set { ch_vcf_final }
    }else{
        ch_vcf_allpos.set { ch_vcf_final }
    }

    FINAL_VCF(
        ch_vcf_final
    )
    
    ch_vcf_final
        .combine(ch_ref_index)
        .dump(tag: 'vcf_and_ref')
        .set { ch_vcf_and_ref }

    CURATE_CONSENSUS(
        ch_vcf_and_ref
    )
    CURATE_CONSENSUS.out.curated_consensus.dump(tag: 'curated_consensus').set { ch_curated }

    emit: 
    ch_curated
}

/*
========================================================================================
    THE END
========================================================================================
*/
