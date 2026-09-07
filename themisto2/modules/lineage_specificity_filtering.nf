process LINEAGE_SPECIFICITY_FILTER {
    tag "${meta.ID}"
    label 'cpu_4'
    label 'mem_8'
    label 'time_queue_from_normal'

    container 'quay.io/sangerpathogens/pandas:2.2.1'

    // Runs once per SPECIES (meta.ID), not per lineage -- one streaming pass over
    publishDir mode: 'copy', path: "${params.outdir}/candidate_marker_filtering/",
               saveAs: { fn -> "${meta.ID}_${fn}" }

    input:
    tuple val(meta), path(unitigs), path(color_sets), path(export_metadata), path(label_mapping)

    output:
    tuple val(meta), path("*_candidate_unitigs.fasta"), emit: unitigs,      optional: true
    tuple val(meta), path("*_stats.txt"),               emit: stats,        optional: true
    tuple val(meta), path("*_specificity.tsv"),         emit: specificity,  optional: true

    script:
    def tg = (meta.target_groups ?: '').trim()
    def lineages_arg = tg ? "--lineages ${tg.tokenize(',').join(' ')}" : "--all-lineages"
    def outside_arg = params.specificity_max_outside != null ? "--max-outside ${params.specificity_max_outside}" : ""
    """
    ${moduleDir}/../bin/lineage_specificity_filter.py \\
        --unitigs ${unitigs} \\
        --color-sets ${color_sets} \\
        --export-metadata ${export_metadata} \\
        --label-mapping ${label_mapping} \\
        ${lineages_arg} \\
        --output-dir . \\
        --stats-output-dir . \\
        --min-freq ${params.candidate_min_freq} \\
        --min-genome-count ${params.candidate_min_genome_count} \\
        --min-lineage-size ${params.candidate_min_genome_count} \\
        --threads ${task.cpus} \\
        ${outside_arg}
    """
}

process CANDIDATE_COLOR_LIST {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_2'
    label 'time_30m'

    container 'quay.io/sangerpathogens/pandas:2.2.1'

    input:
    tuple val(meta), path(candidate_fasta)

    output:
    tuple val(meta), path(file_colors), emit: file_colors

    script:
    file_colors = "candidate_file_colors.txt"
    """
    readlink -f ${candidate_fasta} > ${file_colors}
    """
}
