process THEMISTO2_BUILD {
    tag "${meta.ID}"
    label 'cpu_32'
    label 'time_queue_from_normal'

    // Set directly here, not via a config withName override -- same meta.stage split as
    // GGCAT (see ggcat.nf): THEMISTO2_BUILD_SPECIES measured ~12.3 GB at ~42k genomes,
    // 16 GB baseline leaves headroom; THEMISTO2_BUILD_CANDIDATE's index is tiny.
    memory = { meta.stage == 'candidate' ? (8.GB * task.attempt) : (16.GB * task.attempt) }

    // request /tmp only if /tmp is actually the temp dir (assumes TMPDIR unset)
    if (!params.temp_dir || params.temp_dir.startsWith("/tmp")) {
        label 'request_temp'
    }
    scratch true

    container "quay.io/sangerpathogens/themisto2:0.0.1"

    publishDir mode: 'copy', path: "${params.outdir}/themisto2/${meta.stage ? "${meta.stage}_" : ''}${meta.ID}_build/"

    input:
    tuple val(meta), path(file_colors_input), path(sbwt_index), path(lcs_index)

    output:
    tuple val(meta), path(index_thm2), emit: index

    script:
    index_thm2 = "index.thm2"
    index_build_params = "--file-colors ${file_colors_input} -o ${index_thm2} -s ${sbwt_index} -l ${lcs_index} -k ${params.color_index_kmer_size} -t ${task.cpus}"

    if (params.temp_dir) {
        temp_storage_location = "${params.temp_dir}/themisto2/${meta.ID}"
        index_build_params += " --temp-dir ${temp_storage_location}"
    } else {
        temp_storage_location = "\$PWD"
        index_build_params += " --temp-dir ${temp_storage_location}"
    }


    """
    mkdir -p ${temp_storage_location}
    sed -i '/^\s*\$/d' "${file_colors_input}"    # Remove blank lines
    themisto2 build ${index_build_params}
    """
}

process THEMISTO2_STATS {
    tag "${meta.ID}"
    label 'cpu_4'
    label 'time_30m'

    // real peak RSS <35 MB at candidate scale; species scale needs real headroom over a
    // ~42k-genome index. Same meta.stage split as THEMISTO2_BUILD above.
    memory = { meta.stage == 'candidate' ? (2.GB * task.attempt) : (16.GB * task.attempt) }

    container "quay.io/sangerpathogens/themisto2:0.0.1"

    input:
    tuple val(meta), path(index_thm2)

    output:
    tuple val(meta), path(index_thm2), emit: index

    script:
    """
    themisto2 stats -i ${index_thm2} -t ${task.cpus}
    """
}

process THEMISTO2_EXPORT {
    tag "${meta.ID}"
    label 'cpu_4'
    label 'time_queue_from_normal'

    // Species-wide only -- no candidate alias exists (candidate index is QC/dump only,
    // never exported), so no meta.stage split needed here.
    memory = { 16.GB * task.attempt }

    container "quay.io/sangerpathogens/themisto2:0.0.1"

    publishDir mode: 'copy', path: "${params.outdir}/themisto2/${meta.stage ? "${meta.stage}_" : ''}${meta.ID}_export/"

    input:
    tuple val(meta), path(index_thm2)

    output:
    tuple val(meta), path("export.unitigs.fa"),    emit: unitigs
    tuple val(meta), path("export.color_sets.txt"), emit: color_sets
    tuple val(meta), path("export.metadata.txt"),  emit: metadata

    script:
    """
    themisto2 export -i ${index_thm2} -o export -t ${task.cpus}
    """
}

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
    label 'time_queue_from_normal'

    // Set directly here, not via a config withName override: ATB-species.thm2 is ~293 GB
    // on disk, far bigger than any generic mem_* label covers -- must route to hugemem.
    // No measured peak yet (unlike the sizes above); revisit once a real run has one.
    cpus   = 8
    memory = { 350.GB * task.attempt }
    queue  = 'hugemem'

    container "quay.io/sangerpathogens/themisto2:0.0.1"

    publishDir mode: 'copy', path: "${params.outdir}/atb_cross_species/${meta.ID}/", enabled: params.publish_intermediate

    input:
    // candidate_fasta: this species/lineage's candidate marker unitigs (the same FASTA
    // FILTER_ATB_MARKERS scores verdicts against downstream).
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
