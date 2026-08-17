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
    // Index A -- always built, species-wide.
    tuple val(meta), path("index_species/file_colors_input.txt"), emit: file_colors
    tuple val(meta), path("index_species/label_mapping.tsv"),      emit: label_mapping
    tuple val(meta), path("index_species/stats.json"),             emit: stats
    // Index B -- one subdirectory per requested lineage, only present when
    // target_groups is set. Not yet consumed downstream (GGCAT/SBWT/Themisto2
    // still only build from Index A) -- wiring that up is separate work.
    tuple val(meta), path("index_target_group/*"), emit: target_group_indexes, optional: true

    script:
    // color_mapping.py takes --assembly-dir (+ --assembly-suffix) OR --assembly-paths,
    // never both -- sniff which one assembly_input actually is (directory vs txt file
    // of paths), rather than relying on which CLI param the user set (there's only one now).
    def assembly_arg = assembly_input.isDirectory() \
        ? "--assembly-dir ${assembly_input} --assembly-suffix ${params.assembly_suffix}" \
        : "--assembly-paths ${assembly_input}"
    // Index B: also build lineage-scoped indexes for user-selected groups (e.g.
    // GPSC1,GPSC2), written alongside the species-wide index. Empty by default,
    // which color_mapping.py treats as "species-wide only".
    def target_groups_arg = params.target_groups ? "--target-groups ${params.target_groups}" : ""
    """
    color_mapping.py \\
        --metadata ${metadata} \\
        --sample-col ${params.sample_col} \\
        --label-col ${params.label_col} \\
        ${assembly_arg} \\
        ${target_groups_arg} \\
        --output_dir .
    """
}
