process COLOR_MAPPING {
    tag "${meta.ID}"
    label 'cpu_4'
    label 'mem_4'
    label 'time_queue_from_normal'

    container 'quay.io/sangerpathogens/pandas:2.2.1'

    publishDir mode: 'copy', path: "${params.outdir}/color_mapping/"

    input:
    // label_map: a real map, or assets/NO_LABEL_MAP when the species has none. The
    // placeholder is staged (path inputs can't be empty) but never passed to the script.
    // All four label options are task inputs on purpose: they change which colour each
    // genome gets, so changing one must rebuild the species index.
    tuple val(meta), path(metadata), path(assembly_input), val(label_missing), val(label_multi), path(label_map), val(unclassified_genomes)

    output:
    tuple val(meta), path("${meta.ID}_file_colors_input.txt"), emit: file_colors
    tuple val(meta), path("${meta.ID}_label_mapping.tsv"),     emit: label_mapping
    tuple val(meta), path("${meta.ID}_stats.json"),            emit: stats
    tuple val(meta), path("${meta.ID}_dropped_unclassified.tsv"), emit: dropped_unclassified, optional: true

    script:
    def assembly_arg = assembly_input.isDirectory() \
        ? "--assembly-dir ${assembly_input} --assembly-suffix ${params.assembly_suffix}" \
        : "--assembly-paths ${assembly_input} --assembly-suffix ${params.assembly_suffix}"
    // Single-quoted for the shell; any ' inside becomes '\'' so it can't break the command.
    def label_missing_arg = label_missing ? "--label-missing '${label_missing.replace("'", "'\\''")}'" : ""
    def label_map_arg = label_map.name != 'NO_LABEL_MAP' ? "--label-map ${label_map}" : ""
    """
    ${moduleDir}/../bin/color_mapping.py \\
        --metadata ${metadata} \\
        --species-name ${meta.ID} \\
        --sample-col ${params.sample_col} \\
        --group-label ${params.group_label} \\
        ${assembly_arg} \\
        ${label_missing_arg} \\
        --label-multi ${label_multi} \\
        ${label_map_arg} \\
        --unclassified ${unclassified_genomes} \\
        --output_dir .
    """
}
