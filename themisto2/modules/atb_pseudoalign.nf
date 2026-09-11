// ATB cross-species check (marker_filtering.nf) -- pipeline-workflow replacement for the
// old bg_excl/markers sbwt set-diff (PAT-3570 showed sbwt difference is colour-blind and
// doesn't scale at species-index level). Pseudoaligns a species' candidate marker unitigs
// against Jarno's ATB-species.thm2 (every named bacterial species ATB has assembled,
// ~12,700 colours) and emits one JSONL record per marker with EVERY colour it touches,
// plus (--report-hit-counts) the per-colour k-mer hit count. FILTER_ATB_MARKERS
// (filter_atb_markers.nf) turns that into the actual PASS/FLAG/ABSENT verdict.
//
// --threshold 0.01 --denominator all is a maximally-loose INTAKE gate, not the pass/fail
// rule: it only decides which colours are low enough to bother reporting at all, so
// nothing anywhere near FILTER_ATB_MARKERS' pass/fail line gets silently dropped before
// it's ever seen. The real >=min-within / <=max-outside decision happens downstream.
process THEMISTO2_ATB_PSEUDOALIGN {
    tag "${meta.ID}"
    label 'cpu_8'
    label 'mem_8'   // placeholder -- real sizing is a withName override in
                     // marker_filtering.config: ATB-species.thm2 is ~293 GB on disk
                     // and must route to hugemem, far bigger than any mem_* label covers.
    label 'time_queue_from_normal'

    container "quay.io/sangerpathogens/themisto2:0.0.1"

    publishDir mode: 'copy', path: "${params.outdir}/atb_cross_species/${meta.ID}/", enabled: params.publish_intermediate

    input:
    // candidate_fasta: this species/lineage's candidate marker unitigs (the same FASTA
    // FILTER_ATB_MARKERS scores verdicts against below).
    tuple val(meta), path(candidate_fasta)
    // atb_index: the ~293 GB ATB-species.thm2 (params.atb_index), staged as a real input
    // (not interpolated from params directly) so container binding/caching is Nextflow's
    // problem, not an assumed singularity bind -- same convention the old bg_index used.
    // Same file for every task: broadcast it as a value channel at the call site.
    path(atb_index)

    output:
    tuple val(meta), path(jsonl), emit: jsonl

    script:
    jsonl = "${meta.ID}_atb_pseudoalign.jsonl"
    """
    echo "${candidate_fasta}" > querylist.txt
    echo "${jsonl}"           > outlist.txt

    themisto2 threshold-pseudoalign -i ${atb_index} \\
        --query-list querylist.txt \\
        --query-output-list outlist.txt \\
        --threshold 0.01 --denominator all --min-hits 1 --report-hit-counts \\
        -t ${task.cpus}
    """
}
