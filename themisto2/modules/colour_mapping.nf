include { species_outdir } from './publish_paths.nf'

process COLOUR_MAPPING {
    tag "${meta.ID}"
    label 'cpu_4'
    label 'mem_4'
    label 'time_queue_from_normal'

    container 'quay.io/sangerpathogens/pandas:2.2.1'

    // Final: the ordered assembly list and the Sample_ID -> group mapping. The QC files
    // (stats.json, dropped_unclassified.tsv) only with --publish_intermediate.
    publishDir mode: 'copy', path: "${species_outdir(meta)}/colour_mapping",
               pattern: "*_{file_colours_input.txt,label_mapping.tsv}"
    publishDir mode: 'copy', path: "${species_outdir(meta)}/colour_mapping",
               pattern: "*_{stats.json,dropped_unclassified.tsv}", enabled: params.publish_intermediate

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
