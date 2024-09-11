include { BAKTA               } from '../modules/bakta.nf'
include { COMBINE_ANNOTATIONS } from '../modules/combine_annotations.nf'

workflow ANNOTATE_BAKTA {
    take:
    assembly_channel //tuple val(meta), path(assembly)
    pre_generated_annotation_channel

    main:
    
    BAKTA(assembly_channel)
    
    if (params.combine_annotations) {
        BAKTA.out.gff.join(pre_generated_annotation_channel) //acts as a self filter if there are no annotations given nothing to join and doesn't emit anything
        | COMBINE_ANNOTATIONS
    }

    emit:
    BAKTA.out.gff
}