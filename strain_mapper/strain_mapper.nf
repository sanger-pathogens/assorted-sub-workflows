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
include { CONVERT_TO_BAM; SAMTOOLS_SORT; INDEX_REF; SAMTOOLS_INDEX_BAM } from './modules/samtools'
include {
    BCFTOOLS_CALL;
    BCFTOOLS_MPILEUP;
    BCFTOOLS_FILTERING;
    BCFTOOLS_RAW_VCF;
    BCFTOOLS_FINAL_VCF;
} from './modules/bcftools'
include { PICARD_MARKDUP; PICARD_ADD_READGROUP } from './modules/picard'
include {
    GATK_REF_DICT;
    GATK_HAPLOTYPECALLER;
    GATK_FILTERING;
    GATK_RAW_VCF;
    GATK_FINAL_VCF;
} from './modules/gatk'
include { CURATE_CONSENSUS } from './modules/curate'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow GATK_WORKFLOW {
    take:
    dedup_reads
    ch_ref_index

    main:
    ch_ref_index
        .map { it -> it[0] }
        .set { reference }

    GATK_REF_DICT(
        reference
    )

    ch_ref_index
        .combine(GATK_REF_DICT.out.ref_dict)
        .dump(tag: 'ch_ref_dict')
        .set { ch_ref_dict }

    // Add artificial read group to satisfy haplotypecaller
    PICARD_ADD_READGROUP(
        dedup_reads
    )

    SAMTOOLS_INDEX_BAM(
        PICARD_ADD_READGROUP.out.rg_added_reads
    )
    SAMTOOLS_INDEX_BAM.out.bam_index
        .combine(ch_ref_dict)
        .dump(tag: 'haplotypecaller_input')
        .set { haplotypecaller_input }

    GATK_HAPLOTYPECALLER(
        haplotypecaller_input
    )
    GATK_HAPLOTYPECALLER.out.vcf
        .dump(tag: 'gatk_vcf_allpos')
        .set { ch_vcf_allpos }

    if (params.keep_raw_vcf){
        GATK_RAW_VCF(
            ch_vcf_allpos
        )
    }

    if (!params.skip_filtering) {
        GATK_FILTERING(
            ch_vcf_allpos
        )
        GATK_FILTERING.out.set { ch_filtered_vcf }
    } else {
        ch_vcf_allpos.set { ch_filtered_vcf }
    }

    emit:
    ch_filtered_vcf
}

workflow BCFTOOLS_WORKFLOW {
    take:
    dedup_reads_and_ref

    main:
    BCFTOOLS_MPILEUP(
        dedup_reads_and_ref
    )
    BCFTOOLS_MPILEUP.out.mpileup_file
        .dump(tag: 'mpileup_file')
        .set { ch_mpileup_file }

    BCFTOOLS_CALL(
        ch_mpileup_file
    )
    BCFTOOLS_CALL.out.vcf_allpos
        .dump(tag: 'vcf_allpos')
        .set { ch_vcf_allpos }

    if (params.keep_raw_vcf && !params.skip_filtering){
        BCFTOOLS_RAW_VCF(
            ch_vcf_allpos
        )
    }

    if (!params.skip_filtering) {
        BCFTOOLS_FILTERING(
            ch_vcf_allpos
        )
        BCFTOOLS_FILTERING.out.set { ch_filtered_vcf }
    } else {
        ch_vcf_allpos.set { ch_filtered_vcf }
    }

    emit:
    ch_filtered_vcf
}

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
        .dump(tag: 'dedup_reads_and_ref')
        .set { dedup_reads_and_ref }

    GATK_WORKFLOW(
        PICARD_MARKDUP.out.dedup_reads,
        ch_ref_index
    )

    GATK_FINAL_VCF(
        GATK_WORKFLOW.out.ch_filtered_vcf
    )

    BCFTOOLS_WORKFLOW(
        dedup_reads_and_ref
    )

    BCFTOOLS_FINAL_VCF(
        BCFTOOLS_WORKFLOW.out.ch_filtered_vcf
    )

    BCFTOOLS_WORKFLOW.out.ch_filtered_vcf
        .combine(ch_ref_index)
        .dump(tag: 'vcf_and_ref')
        .set { ch_vcf_and_ref }

    //TODO We could move consensus creation into BCFTOOLS_WORKFLOW and have that workflow?
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
