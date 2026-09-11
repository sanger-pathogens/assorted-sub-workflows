process CHECKPOINT_COUNT {
    tag "${stage}:${meta.ID}"
    label 'cpu_4'
    label 'mem_4'
    label 'time_queue_from_normal'

    errorStrategy 'ignore'

    container "quay.io/sangerpathogens/themisto2:0.0.1"

    publishDir mode: 'copy', path: "${params.outdir}/checkpoints/", enabled: params.publish_intermediate

    input:
    tuple val(meta), val(order), val(stage), val(kind), path(target)

    output:
    tuple val(meta), path(row_tsv), emit: row

    script:
    row_tsv = "${stage}.${meta.ID}.checkpoint.tsv"
    def species = meta.species ?: meta.ID
    def order_key = String.format('%03d', order as int)
    """
    n_kmers=""; n_colors=""; n_unitigs=""; n_seqs=""; sum_bp=""
    min_len=""; median_len=""; max_len=""; n_revcomp_dupes=""

    case "${kind}" in
      colorfile)
        n_seqs=\$(grep -c . "${target}" || true)
        ;;
      fasta)
        n_seqs=\$(grep -c '^>' "${target}" || true)
        sum_bp=\$(grep -v '^>' "${target}" | tr -d '\\n' | wc -c || true)
        # Per-record lengths (one number per line, blank input -> no output), sorted,
        # then min/median/max off that sorted list -- awk only, no seqkit dependency
        # (the themisto2 container doesn't have it).
        awk '/^>/{if(seqlen){print seqlen}; seqlen=0; next} {seqlen+=length(\$0)} END{if(seqlen){print seqlen}}' \\
            "${target}" | sort -n > lengths.txt
        if [ -s lengths.txt ]; then
            min_len=\$(head -n1 lengths.txt)
            max_len=\$(tail -n1 lengths.txt)
            median_len=\$(awk '{a[NR]=\$1} END{n=NR; if(n%2==1){print a[(n+1)/2]} else {print int((a[n/2]+a[n/2+1])/2)}}' lengths.txt)
        fi
        # Forward/reverse-complement duplicate check (e.g. SBWT_DUMP_UNITIGS reports
        # both strands of a unitig as separate records -- caught this via a manual
        # check on the 7PET candidate markers, PAT-3592). Canonicalise each sequence
        # (seq or its revcomp, whichever sorts first) and count records that collapse
        # onto an already-seen canonical form. O(n) revcomps -- fine at marker/candidate
        # scale (thousands of records) but too slow to run unguarded at species-wide
        # scale (~5M records), so skip above REVCOMP_CHECK_MAX rather than silently
        # eating minutes on every full-species run.
        REVCOMP_CHECK_MAX=200000
        if [ "\${n_seqs:-0}" -gt 0 ] && [ "\${n_seqs}" -le "\${REVCOMP_CHECK_MAX}" ]; then
            n_revcomp_dupes=\$(awk '
                function revcomp(s,    i, c, rc) {
                    rc = ""
                    for (i = length(s); i >= 1; i--) {
                        c = substr(s, i, 1)
                        if (c == "A") rc = rc "T"
                        else if (c == "T") rc = rc "A"
                        else if (c == "C") rc = rc "G"
                        else if (c == "G") rc = rc "C"
                        else rc = rc c
                    }
                    return rc
                }
                /^>/ { if (seq != "") { canon = (seq < revcomp(seq)) ? seq : revcomp(seq); if (canon in seen) dupes++; else seen[canon] = 1 }; seq = ""; next }
                { seq = seq \$0 }
                END { if (seq != "") { canon = (seq < revcomp(seq)) ? seq : revcomp(seq); if (canon in seen) dupes++; else seen[canon] = 1 }; print dupes+0 }
            ' "${target}")
        else
            n_revcomp_dupes="skipped_gt_\${REVCOMP_CHECK_MAX}"
        fi
        ;;
      themisto)
        themisto2 stats -i "${target}" -t ${task.cpus} > stats.txt
        n_kmers=\$(sed -n 's/^Number of k-mers: //p' stats.txt)
        n_colors=\$(sed -n 's/^Number of colors: //p' stats.txt)
        n_unitigs=\$(sed -n 's/^Number of forward unitigs (not bidirected): //p' stats.txt)
        ;;
      tsv)
        # Generic header+rows TSV (e.g. design_primers.py's *_primers.tsv /
        # *_no_primers.tsv) -- n_seqs = data row count, header excluded.
        n_seqs=\$(tail -n +2 "${target}" | grep -c . || true)
        ;;
      *)
        echo "CHECKPOINT_COUNT: unknown kind '${kind}'" >&2; exit 1
        ;;
    esac

    printf 'order\\tstage\\tid\\tspecies\\tn_kmers\\tn_colors\\tn_unitigs\\tn_seqs\\tsum_bp\\tmin_len\\tmedian_len\\tmax_len\\n' > "${row_tsv}"
    printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \\
        "${order_key}" "${stage}" "${meta.ID}" "${species}" \\
        "\${n_kmers}" "\${n_colors}" "\${n_unitigs}" "\${n_seqs}" "\${sum_bp}" \\
        "\${min_len}" "\${median_len}" "\${max_len}" >> "${row_tsv}"
    """
}

process CHECKPOINT_COUNT_SBWT {
    tag "${stage}:${meta.ID}"
    label 'cpu_4'
    label 'mem_2'
    label 'time_30m'

    errorStrategy 'ignore'

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
