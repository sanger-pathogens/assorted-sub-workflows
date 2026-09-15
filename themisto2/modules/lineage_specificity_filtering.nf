process LINEAGE_SPECIFICITY_FILTER {
    tag "${meta.ID}"
    label 'cpu_4'
    label 'mem_8'
    label 'time_queue_from_normal'

    container 'quay.io/sangerpathogens/pandas:2.2.1'

    // Runs once per SPECIES (meta.ID), not per lineage -- one streaming pass over
    // the species-wide export scores every --target_groups lineage at once and
    // emits one candidate FASTA per lineage (globbed below).
    publishDir mode: 'copy', path: "${params.outdir}/candidate_filter/${meta.ID}/"

    input:
    // The SPECIES-wide Themisto2 export + species label_mapping (all genomes,
    // every lineage) -- not any one lineage's own export. Needed so the filter can
    // see outside-lineage presence, not just within-lineage.
    tuple val(meta), path(unitigs), path(color_sets), path(export_metadata), path(label_mapping)

    output:
    // Data-driven filenames ('{lineage_id}_...') -- glob, one set per requested lineage.
    tuple val(meta), path("*_candidate_unitigs.fasta"), emit: unitigs
    tuple val(meta), path("*_stats.txt"),               emit: stats
    tuple val(meta), path("*_specificity.tsv"),         emit: specificity

    script:
    def lineages = meta.target_groups.tokenize(',').join(' ')
    def outside_arg = params.specificity_max_outside != null ? "--max-outside ${params.specificity_max_outside}" : ""
    """
    ${moduleDir}/../bin/lineage_specificity_filter.py \\
        --unitigs ${unitigs} \\
        --color-sets ${color_sets} \\
        --export-metadata ${export_metadata} \\
        --label-mapping ${label_mapping} \\
        --lineages ${lineages} \\
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
