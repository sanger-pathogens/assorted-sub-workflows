// Pass/fail half of the ATB cross-species check (see atb_pseudoalign.nf). Reads the JSONL
// THEMISTO2_ATB_PSEUDOALIGN produced and, per marker, scores hit_frac[colour] =
// kmer_hits[colour] / marker_kmer_count against every ATB species colour: PASS iff the
// target species' hit_frac >= atb_min_within AND every other (non-excluded) colour's
// hit_frac <= atb_max_outside. See bin/filter_atb_markers.py's docstring for the full rule,
// including how ATB's lettered species-splits and known-ambiguous colour names are handled.
process FILTER_ATB_MARKERS {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_2'
    label 'time_30m'

    container 'quay.io/sangerpathogens/pandas:2.2.1'

    publishDir mode: 'copy', path: "${params.outdir}/atb_cross_species/${meta.ID}/"

    input:
    // atb_target_species: this species' ATB colour name(s) (manifest atb_target_species
    // column, space-separated) -- NOT part of meta, see manifest_parse.nf.
    tuple val(meta), path(jsonl), path(candidate_fasta), val(atb_target_species)
    // color_names: ATB's color_names.txt (params.atb_color_names), staged as a real input
    // -- same file for every task, broadcast as a value channel at the call site.
    path(color_names)

    output:
    // PASS/FLAG/ABSENT fasta+ids are always written by the script (possibly empty --
    // "no marker passed" is a real, visible result, not a missing file).
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
    ${moduleDir}/../bin/filter_atb_markers.py \\
        --jsonl ${jsonl} \\
        --fasta ${candidate_fasta} \\
        --color-names ${color_names} \\
        --target-species ${atb_target_species} \\
        --min-within ${params.atb_min_within} \\
        --max-outside ${params.atb_max_outside} \\
        --exclude-species ${params.atb_exclude_species.join(' ')} \\
        --out ${prefix}
    """
}
