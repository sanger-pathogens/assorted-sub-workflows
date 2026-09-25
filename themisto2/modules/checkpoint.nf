// Every checkpoint stage, in pipeline order, with the step (process) whose output it counts.
// pipeline_counts.tsv lists rows in this order. To add a stage, insert it here where it
// happens in the pipeline; no numbering to update.
def checkpoint_steps() {
    return [
        species_colourfile           : 'BUILD_COLOUR_INDEX:COLOUR_MAPPING',
        species_ggcat_unitigs        : 'BUILD_COLOUR_INDEX:GGCAT_SPECIES',
        species_sbwt_unitigs         : 'BUILD_COLOUR_INDEX:SBWT_DUMP_UNITIGS_SPECIES',
        species_themisto_index       : 'BUILD_COLOUR_INDEX:THEMISTO2_BUILD_SPECIES',
        species_export_unitigs       : 'BUILD_COLOUR_INDEX:THEMISTO2_EXPORT_SPECIES',
        candidate_specificity_filter : 'MARKER_FILTERING:LINEAGE_SPECIFICITY_FILTER',
        candidate_ggcat_unitigs      : 'MARKER_FILTERING:GGCAT_CANDIDATE',
        candidate_sbwt_unitigs       : 'MARKER_FILTERING:SBWT_DUMP_UNITIGS_CANDIDATE',
        candidate_themisto_index     : 'MARKER_FILTERING:THEMISTO2_BUILD_CANDIDATE',
        candidate_export_unitigs     : 'MARKER_FILTERING:THEMISTO2_EXPORT_CANDIDATE',
        markers_atb_checked          : 'MARKER_FILTERING:FILTER_ATB_MARKERS',
        markers_postproc_pass        : 'POST_PROCESS_MARKERS',
        markers_postproc_reject      : 'POST_PROCESS_MARKERS',
    ]
}

def CHECKPOINT_HEADER = 'step\\tstage\\tid\\tspecies\\tn_kmers\\tn_colours\\tn_unitigs\\tn_seqs\\tsum_bp\\tmin_len\\tmedian_len\\tmax_len\\tn_revcomp_dupes'
def CHECKPOINT_ROW_FMT = '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s'

// The step (process) a stage counts, e.g. BUILD_COLOUR_INDEX:GGCAT_SPECIES.
def checkpoint_step(key) {
    def found = checkpoint_steps().get(key)
    if (found == null) {
        throw new IllegalArgumentException("Unknown checkpoint stage '${key}': add it to checkpoint_steps() in checkpoint.nf")
    }
    return found
}

// Row file name <position>.<stage>.<id>.checkpoint.tsv. The position in checkpoint_steps()
// only orders the rows in CHECKPOINT_REPORT; it never appears in the report.
def checkpoint_row_name(key, id) {
    checkpoint_step(key)
    def position = checkpoint_steps().keySet().toList().indexOf(key)
    return "${String.format('%02d', position)}.${key}.${id}.checkpoint.tsv"
}

process CHECKPOINT_FASTA {
    tag "${stage_name}:${meta.ID}"
    label 'cpu_4'
    label 'mem_4'
    label 'time_queue_from_normal'

    // No errorStrategy here: nextflow-commons' default retries out-of-memory/time kills
    // (exit 130/140) with doubled memory, and ignores any other failure, so a checkpoint
    // never stops the run.

    container 'quay.io/biocontainers/seqkit:2.10.0--h9ee0642_0'

    input:
    tuple val(meta), val(stage_name), val(kind), path(target)

    output:
    tuple val(meta), path(row_tsv), emit: row

    script:
    def step_name = checkpoint_step(stage_name)
    def species = meta.species ?: meta.ID
    row_tsv = checkpoint_row_name(stage_name, meta.ID)
    def k = params.colour_index_kmer_size
    """
    case "${kind}" in
      colourfile)
        n_seqs=\$(grep -c . "${target}" || true)
        ;;

      fasta)
        seqkit stats -a -T "${target}" > stats.tsv
        n_seqs=\$(awk 'NR==2 {print \$4}' stats.tsv)
        sum_bp=\$(awk 'NR==2 {print \$5}' stats.tsv)
        min_len=\$(awk 'NR==2 {print \$6}' stats.tsv)
        median_len=\$(awk 'NR==2 {print \$10}' stats.tsv)
        max_len=\$(awk 'NR==2 {print \$8}' stats.tsv)

        # k-mer positions: sum of (length - k + 1) over records at least k long.
        n_kmers=\$(seqkit fx2tab -n -i -l "${target}" | awk -F'\\t' -v k=${k} '\$2 >= k {s += \$2 - k + 1} END {print s + 0}')

        seqkit rmdup -s -o deduped.fasta -d dupes.fasta "${target}" >/dev/null 2>&1 || true
        n_revcomp_dupes=\$(seqkit stats -T dupes.fasta 2>/dev/null | awk 'NR==2 {print \$4}' || echo 0)
        ;;

      *)
        echo "CHECKPOINT_FASTA: unexpected kind '${kind}' (expected fasta or colourfile)" >&2
        exit 1
        ;;
    esac

    printf '${CHECKPOINT_HEADER}\\n' > "${row_tsv}"
    printf '${CHECKPOINT_ROW_FMT}\\n' \\
        "${step_name}" "${stage_name}" "${meta.ID}" "${species}" \\
        "\${n_kmers:-}" "\${n_colours:-}" "\${n_unitigs:-}" "\${n_seqs:-}" "\${sum_bp:-}" \\
        "\${min_len:-}" "\${median_len:-}" "\${max_len:-}" "\${n_revcomp_dupes:-}" >> "${row_tsv}"
    """
}

process CHECKPOINT_THEMISTO {
    tag "${stage_name}:${meta.ID}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_30m'

    // Same retry behaviour as CHECKPOINT_FASTA (nextflow-commons default).

    container 'quay.io/biocontainers/seqkit:2.10.0--h9ee0642_0'

    // Reads THEMISTO2_STATS' stats.txt rather than loading the index again: loading the
    // species index here ran out of memory (exit 130) on the 7PET run of 23 Sep 2026.
    input:
    tuple val(meta), val(stage_name), path(stats_txt)

    output:
    tuple val(meta), path(row_tsv), emit: row

    script:
    def step_name = checkpoint_step(stage_name)
    def species = meta.species ?: meta.ID
    row_tsv = checkpoint_row_name(stage_name, meta.ID)
    """
    n_kmers=\$(sed -n 's/^Number of k-mers: //p' "${stats_txt}")
    n_colours=\$(sed -n 's/^Number of colors: //p' "${stats_txt}")
    n_unitigs=\$(sed -n 's/^Number of forward unitigs (not bidirected): //p' "${stats_txt}")

    printf '${CHECKPOINT_HEADER}\\n' > "${row_tsv}"
    printf '${CHECKPOINT_ROW_FMT}\\n' \\
        "${step_name}" "${stage_name}" "${meta.ID}" "${species}" \\
        "\${n_kmers:-}" "\${n_colours:-}" "\${n_unitigs:-}" "\${n_seqs:-}" "\${sum_bp:-}" \\
        "\${min_len:-}" "\${median_len:-}" "\${max_len:-}" "\${n_revcomp_dupes:-}" >> "${row_tsv}"
    """
}

// One report per species, rows in checkpoint_steps() order (row files sort by their
// position prefix), always published to results/<species>/checkpoint/.
process CHECKPOINT_REPORT {
    tag "${species}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_30m'

    container 'quay.io/biocontainers/seqkit:2.10.0--h9ee0642_0'

    publishDir mode: 'copy', path: "${params.outdir}/${species}/checkpoint"

    input:
    tuple val(species), path(rows, stageAs: 'rows/*')

    output:
    tuple val(species), path('pipeline_counts.tsv'), emit: report

    script:
    """
    first=1
    for f in \$(ls rows | LC_ALL=C sort); do
        if [ "\$first" = 1 ]; then
            cat "rows/\$f"
            first=0
        else
            tail -n +2 "rows/\$f"
        fi
    done > pipeline_counts.tsv
    """
}
