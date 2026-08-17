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
    tuple val(meta), path("file_colors_input.txt"), emit: file_colors
    tuple val(meta), path("label_mapping.tsv"),      emit: label_mapping
    tuple val(meta), path("stats.txt"),              emit: stats

    script:
    // color_mapping.py takes --assembly-dir (+ --assembly-suffix) OR --assembly-paths,
    // never both -- sniff which one assembly_input actually is (directory vs txt file
    // of paths), rather than relying on which CLI param the user set (there's only one now).
    def assembly_arg = assembly_input.isDirectory() \
        ? "--assembly-dir ${assembly_input} --assembly-suffix ${params.assembly_suffix}" \
        : "--assembly-paths ${assembly_input}"
    // Index B: restrict to user-selected lineages (e.g. GPSC1,GPSC2). Empty by
    // default, which color_mapping.py treats as "keep everything" (Index A).
    def target_groups_arg = params.target_groups ? "--target-groups ${params.target_groups}" : ""
    """
    color_mapping.py \\
        --metadata ${metadata} \\
        --sample-col ${params.sample_col} \\
        --label-col ${params.label_col} \\
        ${assembly_arg} \\
        ${target_groups_arg} \\
        --out-dir .
    """
}
