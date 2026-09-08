

process CHECKPOINT_COUNT {
    tag "${stage}:${meta.ID}"
    label 'cpu_4'
    label 'mem_4'
    label 'time_queue_from_normal'

    // Diagnostic side-channel: a broken counter must never fail a real run.
    errorStrategy 'ignore'

    container "quay.io/sangerpathogens/themisto2:0.0.1"

    publishDir mode: 'copy', path: "${params.outdir}/checkpoints/", enabled: params.publish_intermediate

    input:
    // order = ascending sort key (see file header); kind = 'colorfile' (one line
    // per genome) | 'fasta' | 'themisto' (.thm2)
    tuple val(meta), val(order), val(stage), val(kind), path(target)

    output:
    tuple val(meta), path(row_tsv), emit: row

    script:
    row_tsv = "${stage}.${meta.ID}.checkpoint.tsv"
    def species = meta.species ?: meta.ID
    def order_key = String.format('%03d', order as int)
    """
    n_kmers=""; n_colors=""; n_unitigs=""; n_seqs=""; sum_bp=""

    case "${kind}" in
      colorfile)
        n_seqs=\$(grep -c . "${target}" || true)
        ;;
      fasta)
        n_seqs=\$(grep -c '^>' "${target}" || true)
        sum_bp=\$(grep -v '^>' "${target}" | tr -d '\\n' | wc -c || true)
        ;;
      themisto)
        themisto2 stats -i "${target}" -t ${task.cpus} > stats.txt
        n_kmers=\$(sed -n 's/^Number of k-mers: //p' stats.txt)
        n_colors=\$(sed -n 's/^Number of colors: //p' stats.txt)
        n_unitigs=\$(sed -n 's/^Number of forward unitigs (not bidirected): //p' stats.txt)
        ;;
      *)
        echo "CHECKPOINT_COUNT: unknown kind '${kind}'" >&2; exit 1
        ;;
    esac

    printf 'order\\tstage\\tid\\tspecies\\tn_kmers\\tn_colors\\tn_unitigs\\tn_seqs\\tsum_bp\\n' > "${row_tsv}"
    printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \\
        "${order_key}" "${stage}" "${meta.ID}" "${species}" \\
        "\${n_kmers}" "\${n_colors}" "\${n_unitigs}" "\${n_seqs}" "\${sum_bp}" >> "${row_tsv}"
    """
}

process CHECKPOINT_COUNT_SBWT {
    tag "${stage}:${meta.ID}"
    label 'cpu_4'
    label 'mem_2'
    label 'time_30m'

    // Diagnostic side-channel: a broken counter must never fail a real run.
    errorStrategy 'ignore'

    // Same bug-fixed local .sif as the rest of the SBWT_* processes (see sbwt.nf).
    container "/data/pam/installs/packages/sbwt-rs-cli/bug_fix_setdiff_commit_f93d92_2026.08.04.13.38.59/sbwt-rs-cli-0.4.2-f93d92c/image/sbwt-rs-cli_bug_fix_setdiff_commit_f93d92_2026.08.04.13.38.59.sif"

    publishDir mode: 'copy', path: "${params.outdir}/checkpoints/rows/", enabled: params.publish_intermediate

    input:
    tuple val(meta), val(order), val(stage), path(sbwt_index)

    output:
    tuple val(meta), path(row_tsv), emit: row

    script:
    row_tsv = "${stage}.${meta.ID}.checkpoint.tsv"
    def species = meta.species ?: meta.ID
    def order_key = String.format('%03d', order as int)
    // 'sbwt check' prints e.g. "Index loaded: 7610122 sets, 7578682 k-mers, k=31"
    // to stderr; parse the k-mer count from there.
    """
    n_kmers=""; n_colors=""; n_unitigs=""; n_seqs=""; sum_bp=""

    sbwt check -i "${sbwt_index}" -t ${task.cpus} > check.log 2>&1 || true
    n_kmers=\$(grep -oE '[0-9]+ k-mers' check.log | grep -oE '[0-9]+' | head -1 || true)

    printf 'order\\tstage\\tid\\tspecies\\tn_kmers\\tn_colors\\tn_unitigs\\tn_seqs\\tsum_bp\\n' > "${row_tsv}"
    printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \\
        "${order_key}" "${stage}" "${meta.ID}" "${species}" \\
        "\${n_kmers}" "\${n_colors}" "\${n_unitigs}" "\${n_seqs}" "\${sum_bp}" >> "${row_tsv}"
    """
}
