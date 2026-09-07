#!/usr/bin/env nextflow

include { BOWTIE2_INDEX } from '../modules/bowtie2'
include { BWA_INDEX } from '../modules/bwa'
include { INDEX_REF as SAMTOOLS_INDEX_REF } from '../modules/samtools'


workflow INDEX_REF { 
    take:
    reference

    main:

    if (params.mapper == "bowtie2") {

        // BOWTIE2 INDEX: reuse a complete index sitting next to the reference,
        // build one for every other reference.
        // A complete index is the 6 files written by bowtie2-build, named on the
        // same prefix BOWTIE2_INDEX uses (reference.baseName).
        reference
        .map { refe -> [ refe, files("${refe.parent}/${refe.baseName}{.1,.2,.3,.4,.rev.1,.rev.2}.bt2*") ] }
        .branch { refe, bt2_index_files ->
                has_index:   bt2_index_files.size() == 6
                needs_index: true
        }
        .set { ch_bt2 }

        BOWTIE2_INDEX( ch_bt2.needs_index.map { refe, bt2_index_files -> refe } )
        .bt2_index
        .mix( ch_bt2.has_index )
        .set { ch_bt2_index }

        ch_bwa_index = Channel.empty()

    } else if (params.mapper == "bwa") {

        // BWA INDEX: reuse a complete index sitting next to the reference,
        // build one for every other reference
        reference
        .map { refe -> [ refe, files("${refe}{.amb,.ann,.bwt,.pac,.sa}") ] }
        .branch { refe, bwa_index_files ->
                has_index:   bwa_index_files.size() == 5
                needs_index: true
        }
        .set { ch_bwa }

        BWA_INDEX( ch_bwa.needs_index.map { refe, bwa_index_files -> refe } )
        .bwa_index
        .mix( ch_bwa.has_index )
        .set { ch_bwa_index }

        ch_bt2_index = Channel.empty()

    } else {
        error "supplied mapper: ${params.mapper} is not currently supported"
    }

    // SAMTOOLS INDEX REF FASTA FOR DOWNSTREAM PROCESSES
    reference
    .branch{ refe ->
            has_faidx: file("${refe}.fai").isFile()
            needs_faidx:  true
    }
    .set { ch_ref }
    
    ch_ref.has_faidx
    .map { refe -> [ refe, file("${refe}.fai") ] }
    | set { ch_ref_preindex }

    SAMTOOLS_INDEX_REF(ch_ref.needs_faidx)
    .mix(ch_ref_preindex)
    | set { ch_ref_index }

    emit:
    ch_bt2_index
    ch_bwa_index
    ch_ref_index
}