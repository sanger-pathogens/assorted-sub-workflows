process FILTER_ATB_MARKERS {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_2'
    label 'time_30m'

    container 'quay.io/sangerpathogens/pandas:2.2.1'

    publishDir mode: 'copy', path: "${params.outdir}/atb_cross_species/${meta.ID}/"

    input:
    tuple val(meta), path(jsonl), path(candidate_fasta), val(atb_target_species), val(atb_exclude_species)
    path(colour_names)

    output:
    tuple val(meta), path("${prefix}_PASS.fasta"),         emit: pass
    tuple val(meta), path("${prefix}_FLAG.fasta"),         emit: flag
    tuple val(meta), path("${prefix}_ABSENT.fasta"),       emit: absent
    tuple val(meta), path("${prefix}_validation.tsv"),     emit: validation
    tuple val(meta), path("${prefix}_species_detail.tsv"), emit: species_detail, optional: true
    tuple val(meta), path("${prefix}_warnings.tsv"),       emit: warnings,       optional: true
    tuple val(meta), path("${prefix}_summary.txt"),        emit: summary

    script:
    prefix = "${meta.ID}_atb_check"
    """
    ${moduleDir}/../bin/atb_cross_species_filter.py \\
        --jsonl ${jsonl} \\
        --fasta ${candidate_fasta} \\
        --color-names ${colour_names} \\
        --target-species ${atb_target_species} \\
        --min-within ${params.atb_min_within} \\
        --max-outside ${params.atb_max_outside} \\
        --exclude-species ${atb_exclude_species} \\
        --out ${prefix}
    """
}
