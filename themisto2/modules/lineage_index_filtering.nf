process CANDIDATE_FILTER {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_2'
    label 'time_30m'

    container 'quay.io/sangerpathogens/pandas:2.2.1'

    publishDir mode: 'copy', path: "${params.outdir}/candidate_filter/${meta.ID}/"

    input:
    // unitigs/color_sets/export_metadata: this lineage's OWN Themisto2 export
    // (BUILD_COLOR_INDEX's lineage branch), not the species-wide one -- every
    // colour in it already is a lineage genome. label_mapping: this lineage's own
    // color_mapping.py output (index_target_group/<group>/label_mapping.tsv), used only
    // to read back the lineage ID.
    tuple val(meta), path(unitigs), path(color_sets), path(export_metadata), path(label_mapping)

    output:
    // Filenames are data-driven ('{lineage_id}_{min_freq_label}_...') -- glob rather
    // than name them explicitly, since min_freq_label depends on params.candidate_min_freq
    // (a preset name or a custom 'freqXpY' token) and this process only ever produces
    // one of each per task.
    tuple val(meta), path("*_candidate_unitigs.fasta"), emit: unitigs
    tuple val(meta), path("*_stats.txt"),                emit: stats

    script:
    // Count and fraction filtering. Drop if min_freq is less than N or if the
    // fraction of unitigs shared across the genomes in that lineage group is less
    // than the selected mode chosen.
    //
    // --stats-output-dir is always passed for now -- TODO: make this optional via
    // its own config flag once there's a reason not to always want it.
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
    // GGCAT/THEMISTO2_BUILD expect a colour-list file (genome paths, one per line,
    // each becoming its own colour), not a raw sequences FASTA -- CANDIDATE_FILTER's
    // output is the latter, so wrap its own resolved path as a one-line list,
    // matching color_mapping.py's canonical-absolute-path approach for the same
    // cross-task-directory reason. One colour is correct here: colour/genome-level
    // distinction was already consumed by CANDIDATE_FILTER's own size-based
    // thresholding, nothing left to distinguish by colour at this point.
    tuple val(meta), path(candidate_fasta)

    output:
    tuple val(meta), path(file_colors), emit: file_colors

    script:
    file_colors = "candidate_file_colors.txt"
    """
    readlink -f ${candidate_fasta} > ${file_colors}
    """
}
