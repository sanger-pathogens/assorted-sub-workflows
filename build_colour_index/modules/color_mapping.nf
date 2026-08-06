process COLOR_MAPPING {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_2'
    label 'time_queue_from_small'

    container 'quay.io/sangerpathogens/pandas:2.2.1'

    publishDir mode: 'copy', path: "${params.outdir}/colour_mapping/${meta.ID}/"

    input:
    tuple val(meta), path(metadata_csv), path(assembly_input)

    output:
    tuple val(meta), path("file_colors_input.txt"), emit: file_colors
    tuple val(meta), path("label_mapping.tsv"),      emit: label_mapping
    tuple val(meta), path("stats.txt"),              emit: stats

    script:
    // color_mapping.py takes --assembly-dir (+ --assembly-suffix) OR --assembly-paths,
    // never both
    def assembly_arg = params.assembly_paths \
        ? "--assembly-paths ${assembly_input}" \
        : "--assembly-dir ${assembly_input} --assembly-suffix ${params.assembly_suffix}"
    """
    color_mapping.py \\
        --metadata ${metadata_csv} \\
        --sample-col ${params.sample_col} \\
        --label-col ${params.label_col} \\
        ${assembly_arg} \\
        --out-dir .
    """
}
