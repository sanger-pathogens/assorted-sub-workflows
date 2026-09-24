def CHECKPOINT_HEADER = 'order\\tstage\\tid\\tspecies\\tn_kmers\\tn_colours\\tn_unitigs\\tn_seqs\\tsum_bp\\tmin_len\\tmedian_len\\tmax_len\\tn_revcomp_dupes'
def CHECKPOINT_ROW_FMT = '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s'

process CHECKPOINT_FASTA {
    tag "${stage}:${meta.ID}"
    label 'cpu_4'
    label 'mem_4'
    label 'time_queue_from_normal'

    errorStrategy 'ignore'

    container 'quay.io/biocontainers/seqkit:2.10.0--h9ee0642_0'

    input:
    tuple val(meta), val(order), val(stage), val(kind), path(target)

    output:
    tuple val(meta), path(row_tsv), emit: row

    script:
    row_tsv = "${stage}.${meta.ID}.checkpoint.tsv"
    def species = meta.species ?: meta.ID
    def order_key = String.format('%03d', order as int)
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
        "${order_key}" "${stage}" "${meta.ID}" "${species}" \\
        "\${n_kmers:-}" "\${n_colours:-}" "\${n_unitigs:-}" "\${n_seqs:-}" "\${sum_bp:-}" \\
        "\${min_len:-}" "\${median_len:-}" "\${max_len:-}" "\${n_revcomp_dupes:-}" >> "${row_tsv}"
    """
}

process CHECKPOINT_THEMISTO {
    tag "${stage}:${meta.ID}"
    label 'cpu_4'
    label 'mem_4'
    label 'time_queue_from_normal'

    errorStrategy 'ignore'

    container 'quay.io/sangerpathogens/themisto2:0.0.1'

    input:
    tuple val(meta), val(order), val(stage), val(kind), path(target)

    output:
    tuple val(meta), path(row_tsv), emit: row

    script:
    row_tsv = "${stage}.${meta.ID}.checkpoint.tsv"
    def species = meta.species ?: meta.ID
    def order_key = String.format('%03d', order as int)
    """
    themisto2 stats -i "${target}" -t ${task.cpus} > stats.txt
    n_kmers=\$(sed -n 's/^Number of k-mers: //p' stats.txt)
    n_colours=\$(sed -n 's/^Number of colors: //p' stats.txt)
    n_unitigs=\$(sed -n 's/^Number of forward unitigs (not bidirected): //p' stats.txt)

    printf '${CHECKPOINT_HEADER}\\n' > "${row_tsv}"
    printf '${CHECKPOINT_ROW_FMT}\\n' \\
        "${order_key}" "${stage}" "${meta.ID}" "${species}" \\
        "\${n_kmers:-}" "\${n_colours:-}" "\${n_unitigs:-}" "\${n_seqs:-}" "\${sum_bp:-}" \\
        "\${min_len:-}" "\${median_len:-}" "\${max_len:-}" "\${n_revcomp_dupes:-}" >> "${row_tsv}"
    """
}
