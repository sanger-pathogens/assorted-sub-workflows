process COMBINE_ANNOTATIONS {
    tag "${meta.ID}"
    label "cpu_4"
    label "mem_16"
    label "time_1"

    publishDir mode: 'copy', pattern: "${meta.ID}_merged_annotation.gff3", path: "${params.outdir}/combined_gffs"

    container 'quay.io/sangerpathogens/gffutils:0.13'

    input:
    tuple val(meta), path(main_annotation), path(annotation)

    script:
    merge_annotations = "${projectDir}/assorted-sub-workflows/annotate_bakta/bin/merge_annotations.py"
    """
    ${merge_annotations} --main_file ${main_annotation} --additional_files ${annotation} -o ${meta.ID}_merged_annotation.gff3 -d db
    """
}