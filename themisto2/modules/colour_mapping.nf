process COLOUR_MAPPING {
    tag "${meta.ID}"
    label 'cpu_4'
    label 'mem_4'
    label 'time_queue_from_normal'

    container 'quay.io/sangerpathogens/pandas:2.2.1'

    publishDir mode: 'copy', path: "${params.outdir}/colour_mapping/"

    input:
    tuple val(meta), path(metadata), path(assembly_input)

    output:
    tuple val(meta), path("${meta.ID}_file_colours_input.txt"), emit: file_colours
    tuple val(meta), path("${meta.ID}_label_mapping.tsv"),     emit: label_mapping
    tuple val(meta), path("${meta.ID}_stats.json"),            emit: stats
    tuple val(meta), path("${meta.ID}_dropped_unclassified.tsv"), emit: dropped_unclassified, optional: true

    script:
    def assembly_arg = assembly_input.isDirectory() \
        ? "--assembly-dir ${assembly_input} --assembly-suffix ${params.assembly_suffix}" \
        : "--assembly-paths ${assembly_input} --assembly-suffix ${params.assembly_suffix}"
    """
    ${moduleDir}/../bin/colour_mapping.py \\
        --metadata ${metadata} \\
        --species-name ${meta.ID} \\
        --sample-col ${params.sample_col} \\
        --group-label ${params.group_label} \\
        ${assembly_arg} \\
        --output_dir .
    """
}
