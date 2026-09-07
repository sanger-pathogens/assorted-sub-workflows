#!/usr/bin/env nextflow

include { BOWTIE2_INDEX } from '../modules/bowtie2'
include { BWA_INDEX } from '../modules/bwa'
include { INDEX_REF as SAMTOOLS_INDEX_REF } from '../modules/samtools'


workflow INDEX_REF {
    take:
    reference        // channel: path(reference), deduplicated

    main:

    // Key every reference by its original path. Index processes re-emit the reference
    // from their work directory, so the path itself cannot be used to join indexes back
    // to reads; this val key is not staged and so is stable on both sides of the split
    // between references we index here and references that were already indexed.
    reference
    .map { refe -> [ refe.toString(), refe ] }
    .set { ch_ref_keyed }

    if (params.mapper == "bowtie2") {

        // BOWTIE2 INDEX: reuse a complete index sitting next to the reference,
        // build one for every other reference.
        // A complete index is the 6 files written by bowtie2-build, named on the
        // same prefix BOWTIE2_INDEX uses (reference.baseName).
        ch_ref_keyed
        .map { ref_key, refe -> [ ref_key, refe, files("${refe.parent}/${refe.baseName}{.1,.2,.3,.4,.rev.1,.rev.2}.bt2*") ] }
        .branch { ref_key, refe, bt2_index_files ->
                has_index:   bt2_index_files.size() == 6
                needs_index: true
        }
        .set { ch_bt2 }

        BOWTIE2_INDEX( ch_bt2.needs_index.map { ref_key, refe, bt2_index_files -> [ ref_key, refe ] } )
        .bt2_index
        .mix( ch_bt2.has_index )
        .set { ch_bt2_index }

        ch_bwa_index = Channel.empty()

    } else if (params.mapper == "bwa") {

        // BWA INDEX: reuse a complete index sitting next to the reference,
        // build one for every other reference
        ch_ref_keyed
        .map { ref_key, refe -> [ ref_key, refe, files("${refe}{.amb,.ann,.bwt,.pac,.sa}") ] }
        .branch { ref_key, refe, bwa_index_files ->
                has_index:   bwa_index_files.size() == 5
                needs_index: true
        }
        .set { ch_bwa }

        BWA_INDEX( ch_bwa.needs_index.map { ref_key, refe, bwa_index_files -> [ ref_key, refe ] } )
        .bwa_index
        .mix( ch_bwa.has_index )
        .set { ch_bwa_index }

        ch_bt2_index = Channel.empty()

    } else {
        error "supplied mapper: ${params.mapper} is not currently supported"
    }

    // SAMTOOLS INDEX REF FASTA FOR DOWNSTREAM PROCESSES
    ch_ref_keyed
    .map { ref_key, refe -> [ ref_key, refe, file("${refe}.fai") ] }
    .branch { ref_key, refe, faidx ->
            has_faidx:   faidx.isFile()
            needs_faidx: true
    }
    .set { ch_ref }

    SAMTOOLS_INDEX_REF( ch_ref.needs_faidx.map { ref_key, refe, faidx -> [ ref_key, refe ] } )
    .ref_index
    .mix( ch_ref.has_faidx )
    .set { ch_ref_index }

    emit:
    ch_bt2_index    // channel: [ val(ref_key), path(reference), path(bt2_index_files) ]
    ch_bwa_index    // channel: [ val(ref_key), path(reference), path(bwa_index_files) ]
    ch_ref_index    // channel: [ val(ref_key), path(reference), path(faidx) ]
}