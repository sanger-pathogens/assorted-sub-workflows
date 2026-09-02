#!/usr/bin/env nextflow

include { BOWTIE2_INDEX } from '../modules/bowtie2'
include { BWA_INDEX } from '../modules/bwa'
include { INDEX_REF as SAMTOOLS_INDEX_REF } from '../modules/samtools'


workflow INDEX_REF { 
    take:
    reference

    main:

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
        ch_bwa_index = Channel.empty()

    } else if (params.mapper == "bwa") {

        // BWA INDEX
        bwa_index_files = Path("${reference}.amb")
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
        ch_bt2_index = Channel.empty()

    // SAMTOOLS INDEX REF FASTA FOR DOWNSTREAM PROCESSES
    references
    .map { refe -> 
        faidx = Path("${refe}.fai")
        if (faidx.isFile()) {
            [reference, faidx]
        } else {
            [reference, null]
        }
    }
    | set { ch_ref_preindex }

    } else {
        SAMTOOLS_INDEX_REF(ch_ref_preindex)
        | set { ch_ref_index }
    }


    emmit:
    ch_bt2_index
    ch_bwa_index
    ch_ref_index
}