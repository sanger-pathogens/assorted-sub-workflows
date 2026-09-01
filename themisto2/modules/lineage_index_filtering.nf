process CANDIDATE_FILTER {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_2'
    label 'time_30m'

    container 'quay.io/sangerpathogens/pandas:2.2.1'

    publishDir mode: 'copy', path: "${params.outdir}/candidate_filter/${meta.ID}/"

    input:
    // This lineage's OWN Themisto2 export (not the species-wide one); label_mapping
    // is only read back for the lineage ID.
    tuple val(meta), path(unitigs), path(color_sets), path(export_metadata), path(label_mapping)

    output:
    // Data-driven filenames ('{lineage_id}_{min_freq_label}_...') -- glob, one per task.
    tuple val(meta), path("*_candidate_unitigs.fasta"), emit: unitigs
    tuple val(meta), path("*_stats.txt"),                emit: stats

    script:
    """
    ${moduleDir}/../bin/core_catchall_filter.py \\
        --unitigs ${unitigs} \\
        --color-sets ${color_sets} \\
        --export-metadata ${export_metadata} \\
        --label-mapping ${label_mapping} \\
        --output-dir . \\
        --stats-output-dir . \\
        --min-freq ${params.candidate_min_freq} \\
        --min-genome-count ${params.candidate_min_genome_count}
    """
}

process CANDIDATE_COLOR_LIST {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_2'
    label 'time_30m'

    container 'quay.io/sangerpathogens/pandas:2.2.1'

    input:
    // GGCAT/THEMISTO2_BUILD want a colour-list file (paths, one per line), not a raw
    // FASTA. Wrap the candidate FASTA's absolute path as a one-line, single-colour list.
    tuple val(meta), path(candidate_fasta)

    output:
    tuple val(meta), path(file_colors), emit: file_colors

    script:
    file_colors = "candidate_file_colors.txt"
    """
    readlink -f ${candidate_fasta} > ${file_colors}
    """
}
