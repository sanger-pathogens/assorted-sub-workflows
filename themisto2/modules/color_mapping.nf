process COLOR_MAPPING {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_2'
    label 'time_queue_from_small'

    container 'quay.io/sangerpathogens/pandas:2.2.1'

    publishDir mode: 'copy', path: "${params.outdir}/colour_mapping/${meta.ID}/"

    input:
    tuple val(meta), path(metadata), path(assembly_input)

    output:
    tuple val(meta), path("index_species/species_file_colors_input.txt"), emit: file_colors
    tuple val(meta), path("index_species/species_label_mapping.tsv"),     emit: label_mapping
    tuple val(meta), path("index_species/species_stats.json"),            emit: stats
    tuple val(meta), path("index_target_group/*/*_file_colors_input.txt"), emit: target_group_file_colors,   optional: true
    tuple val(meta), path("index_target_group/*/*_label_mapping.tsv"),     emit: target_group_label_mapping, optional: true
    tuple val(meta), path("index_target_group/*/*_stats.json"),            emit: target_group_stats,         optional: true

    script:
    def assembly_arg = assembly_input.isDirectory() \
        ? "--assembly-dir ${assembly_input} --assembly-suffix ${params.assembly_suffix}" \
        : "--assembly-paths ${assembly_input} --assembly-suffix ${params.assembly_suffix}"
    def target_groups_arg = params.target_groups ? "--target-groups ${params.target_groups}" : ""
    """
    color_mapping.py \\
        --metadata ${metadata} \\
        --sample-col ${params.sample_col} \\
        --group-label ${params.group_label} \\
        ${assembly_arg} \\
        ${target_groups_arg} \\
        --output_dir .
    """
}
