process MDMCLEANER_SKIP {
    // Isolates (pure cultures) are not cleaned: pass the genome through unchanged, under the file name
    // MDMCLEANER would have produced, so the post-QC steps treat both genome types identically.
    tag "${meta.ID}"
    label "cpu_1"
    label "mem_1"
    label "time_30m"

    container  'quay.io/biocontainers/seqkit:2.10.0--h9ee0642_0'

    input:
    tuple val(group_key), val(meta), path(fasta)

    output:
    tuple val(group_key), val(meta), path("*_filtered_kept_contigs.${output_fasta_ext}"), emit: results

    script:
    fasta_ext = params.fasta_ext.replaceAll(/^\./, '')
    output_fasta_ext = "${fasta_ext}.gz"
    name = fasta.name.endsWith(".${fasta_ext}") ? fasta.name[0..<(fasta.name.length() - fasta_ext.length() - 1)] : fasta.baseName
    """
    if gzip -t ${fasta} 2>/dev/null; then
        cp ${fasta} ${name}_filtered_kept_contigs.${output_fasta_ext}
    else
        gzip -c ${fasta} > ${name}_filtered_kept_contigs.${output_fasta_ext}
    fi
    """
}

process MDMCLEANER {
    tag "${meta.ID}"
    label "cpu_16"
    label "mem_8"
    label "time_12"

    container  'quay.io/biocontainers/mdmcleaner:0.8.7--pyh7cba7a3_0'

    publishDir enabled: params.publish_mdmcleaner , mode: 'copy', path: "${params.outdir}/mdmcleaner/"

    input:
    tuple val(group_key), val(meta), path(fasta)

    output:
    tuple val(group_key), val(meta), path("*.${output_fasta_ext}"), emit: results

    script:
    fasta_ext = params.fasta_ext.replaceAll(/^\./, '')
    output_fasta_ext = "${fasta_ext}.gz"
    // genome name = input file name without .<fasta_ext>[.gz]; MDMcleaner names its outputs after it
    genome_name = fasta.name.endsWith('.gz') ? fasta.name[0..<-3] : fasta.name
    genome_name = genome_name.endsWith(".${fasta_ext}") ? genome_name[0..<-(fasta_ext.length() + 1)] : genome_name
    """
    # MDMcleaner splits contig names at ':' when it matches rRNA hits back to contigs, so a name such as
    # '3300029572.a:Ga0244925_108424' (JGI/IMG style) raises KeyError and the genome is dropped by the
    # 'ignore' errorStrategy without any warning (305 El Mouali 2026 genomes, Oct 2026). Swap ':' for a
    # placeholder in a copy of the input, then restore it in the cleaned output, so contig names come out
    # exactly as they went in.
    # barrnap (inside MDMcleaner) greps every line of nhmmer's output for /fail|error/i, and that output names the
    # input file, so a genome whose FILE NAME contains "fail" or "error" is reported as an nhmmer failure and
    # dropped (found 2026-10-16 with a test genome called mag_gunc_fail). Run MDMcleaner on a neutrally named
    # copy and rename its outputs back afterwards.
    mkdir -p header_safe
    if gzip -t ${fasta} 2>/dev/null; then
        zcat ${fasta} | sed '/^>/ s/:/__COLON__/g' | gzip > header_safe/qcgenome.${fasta_ext}.gz
        staged=header_safe/qcgenome.${fasta_ext}.gz
    else
        sed '/^>/ s/:/__COLON__/g' ${fasta} > header_safe/qcgenome.${fasta_ext}
        staged=header_safe/qcgenome.${fasta_ext}
    fi

    mdmcleaner set_configs --db_basedir ${params.mdmcleaner_db} --threads ${task.cpus}
    mdmcleaner clean -i \$staged -o mdmcleaner_output -t ${task.cpus} -c mdmcleaner.config --fast_run

    # restore the genome's own name on MDMcleaner's output folder and files
    if [ -d mdmcleaner_output/qcgenome ]; then
        for file in mdmcleaner_output/qcgenome/qcgenome_*; do
            mv "\$file" "mdmcleaner_output/qcgenome/${genome_name}_\${file##*/qcgenome_}"
        done
        mv mdmcleaner_output/qcgenome mdmcleaner_output/${genome_name}
    fi

    for file in mdmcleaner_output/*/*.fasta.gz; do
        zcat "\$file" | sed '/^>/ s/__COLON__/:/g' | gzip > "\$file.tmp"
        mv "\$file.tmp" "\$file"
    done

    # Rename files to match desired extension as mdmcleaner is always fasta
    for file in mdmcleaner_output/*/*.fasta.gz; do
        new_name=\$(echo "\$file" | sed "s|\\.fasta\\.gz\$|.${output_fasta_ext}|")
        mv "\$file" "\$new_name"
    done

    mv mdmcleaner_output/*/*_filtered_kept_contigs.${output_fasta_ext} .
    """
}

process MDMCLEANER_AUDIT {
    executor 'local'

    publishDir mode: 'copy', path: "${params.outdir}/report/"

    input:
    path(sent)
    path(cleaned)

    output:
    path("mdmcleaner_skipped.tsv"), emit: skipped
    stdout emit: n_skipped

    script:
    """
    printf 'batch\\tgenome\\n' > mdmcleaner_skipped.tsv
    sort -u ${cleaned} > cleaned.sorted
    sort -u ${sent} | awk -F'\\t' 'NR==FNR {done[\$0]=1; next} !(\$0 in done)' cleaned.sorted - >> mdmcleaner_skipped.tsv
    printf '%s' \$(( \$(wc -l < mdmcleaner_skipped.tsv) - 1 ))
    """
}
