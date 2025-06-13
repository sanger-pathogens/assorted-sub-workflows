#!/usr/bin/env nextflow

/*
========================================================================================
CONVERTED FOR SUBWORKFLOW PURPOSES
========================================================================================
*/

// This altered abundance estimation includes changes that relate to the input/output of the pipeline as well as changes to the location
// of module files. In the main body of code there is also an included skip for the sourmash gather and subsetting of GTDB which is
// not relevent to RVI's use case for now - maybe in the future this will be required.


/*
========================================================================================
    IMPORT MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULES
//
include { CLEANUP_SORTED_BAM_FILES; CLEANUP_TRIMMED_FASTQ_FILES; CLEANUP_INSTRAIN_OUTPUT } from '../modules/cleanup.nf'
include { MERGE_FASTQS } from '../modules/merge_fastq.nf'
include { SOURMASH_SKETCH; SOURMASH_GATHER } from '../modules/sourmash.nf'
include { SUBSET_GTDB } from '../modules/subset_fasta.nf'
include { BOWTIE_INDEX; BOWTIE2SAMTOOLS; GET_OVERALL_MAPPING_RATE } from '../modules/bowtie.nf'
include { GENERATE_STB; INSTRAIN_PROFILE; INSTRAIN_QUICKPROFILE; GENERATE_INSTRAIN_SUMMARY } from '../modules/instrain.nf'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow ABUNDANCE_ESTIMATION{
    take:
    fastq_path_ch        // tuple( meta, read_1, read_2 )

    main:
    MERGE_FASTQS(fastq_path_ch)

    if (params.sourmash_subset_abundance_estimation) {

        SOURMASH_SKETCH(MERGE_FASTQS.out.merged_fastq)

        SOURMASH_GATHER(SOURMASH_SKETCH.out.sketch)

        SUBSET_GTDB(SOURMASH_GATHER.out.sourmash_genomes)

        SUBSET_GTDB.out.subset_genome.set{ genomes_channel }

        BOWTIE_INDEX(genomes_channel)

        BOWTIE_INDEX.out.bowtie_index.set{ index_channel }

        GENERATE_STB(SOURMASH_GATHER.out.sourmash_genomes)

        GENERATE_STB.out.stb_ch.set{ stb_channel }
    } else {

        genomes_channel = Channel.fromPath("${params.genome_file_abundance_estimation}")
        
        index_channel = Channel.of("${params.precomputed_index_abundance_estimation}")

        stb_channel = Channel.fromPath("${params.stb_file_abundance_estimation}")
    }

    fastq_path_ch.combine(index_channel).set{ bowtie_mapping_ch }
    BOWTIE2SAMTOOLS(bowtie_mapping_ch, params.bowtie2_samtools_threads_abundance_estimation)

    if (params.cleanup_intermediate_files_abundance_estimation) {
        CLEANUP_TRIMMED_FASTQ_FILES(BOWTIE2SAMTOOLS.out.trimmed_fastqs)
    }

    GET_OVERALL_MAPPING_RATE(BOWTIE2SAMTOOLS.out.map_rate_ch.collect())

    if (!params.bowtie2_samtools_only_abundance_estimation) {
        BOWTIE2SAMTOOLS.out.bam_file
        .combine(stb_channel)
        .combine(genomes_channel)
        .set { instrain_profiling_ch }

        if (!params.instrain_quick_profile_abundance_estimation){
            INSTRAIN_PROFILE(instrain_profiling_ch)

            INSTRAIN_PROFILE.out.genome_info_file
            .collect() { it[1] }
            .set { genome_info_files }

            INSTRAIN_PROFILE.out.meta_workdir
            .set { instrain_meta_workdir_ch } 

            GENERATE_INSTRAIN_SUMMARY(genome_info_files)
        }else{
            INSTRAIN_QUICKPROFILE(instrain_profiling_ch)

            INSTRAIN_QUICKPROFILE.out.meta_workdir
            .set { instrain_meta_workdir_ch } 
        }
    }

    if (params.cleanup_intermediate_files_abundance_estimation && params.bowtie2_samtools_only_abundance_estimation) {
        bam_files = BOWTIE2SAMTOOLS.out.bam_file.map { it[1] }
        CLEANUP_SORTED_BAM_FILES(bam_files)
    }
    else if(params.cleanup_intermediate_files_abundance_estimation) {
        CLEANUP_SORTED_BAM_FILES(INSTRAIN.out.sorted_bam)
    }

    if (params.cleanup_intermediate_files_abundance_estimation && !params.bowtie2_samtools_only_abundance_estimation) {
        CLEANUP_INSTRAIN_OUTPUT(instrain_meta_workdir_ch)
    }
}
