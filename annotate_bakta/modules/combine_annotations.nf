process COMBINE_ANNOTATIONS {
    tag "${meta.ID}"
    label "cpu_4"
    label "mem_16"
    label "time_1"

    publishDir mode: 'copy', pattern: "${meta.ID}.gff3", path: "${params.outdir}/combined"

    container 'quay.io/ssd28/gsoc-experimental/bcf_2_pseudosequence:0.0.2'

    input:
    tuple val(meta), path(main_annotation), path('annotations???')

    script:
    """
    merge_annotations.py --main_file ${main_annotation} --additional_files annotations* -o ${meta.ID}_merged_annotation.gff3 -d ${task.workDir}/db
    """
}