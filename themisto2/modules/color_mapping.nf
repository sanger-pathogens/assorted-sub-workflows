process COLOR_MAPPING {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_2'
    label 'time_30m'

    container 'quay.io/sangerpathogens/pandas:2.2.1'

    publishDir mode: 'copy', path: "${params.outdir}/colour_mapping/${meta.ID}/"

    input:
    tuple val(meta), path(metadata), path(assembly_input)

    output:
    // Species-wide only. Per-lineage colour files used to be built here for Index B;
    // Index B is gone (PAT-3570) and step 07 filtering reads the species-wide export,
    // so meta.target_groups no longer reaches this step -- it selects lineages in
    // LINEAGE_SPECIFICITY_FILTER instead.
    tuple val(meta), path("index_species/species_file_colors_input.txt"), emit: file_colors
    tuple val(meta), path("index_species/species_label_mapping.tsv"),     emit: label_mapping
    tuple val(meta), path("index_species/species_stats.json"),            emit: stats

    script:
    def assembly_arg = assembly_input.isDirectory() \
        ? "--assembly-dir ${assembly_input} --assembly-suffix ${params.assembly_suffix}" \
        : "--assembly-paths ${assembly_input} --assembly-suffix ${params.assembly_suffix}"
    """
    ${moduleDir}/../bin/color_mapping.py \\
        --metadata ${metadata} \\
        --sample-col ${params.sample_col} \\
        --group-label ${params.group_label} \\
        ${assembly_arg} \\
        --output_dir .
    """
}
