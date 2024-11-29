process GUBBINS {
    label 'cpu_16'
    label 'mem_1'
    label 'time_12'

    conda 'bioconda::gubbins=3.2.1'
    container 'quay.io/biocontainers/gubbins:3.2.1--py38pl5321h4c6a040_1'

    publishDir "${params.outdir}/gubbins", mode: 'copy', overwrite: true, pattern: "${prefix}.*"

    input:
    path(msa)

    output:
    tuple path(msa), path(recombination_pred), emit: recombination_pred

    script:
    prefix = "gubbins_out"
    recombination_pred = "${prefix}.recombination_predictions.gff"
    recombination_free_tree = "${prefix}.node_labelled.final_tree.tre"
    """
    run_gubbins.py \\
        --prefix ${prefix} \\
        --filter-percentage ${params.gubbins_filter_percentage} \\
        --threads ${task.cpus} \\
        ${msa}
    """
}

process GUBBINS_MASK {
    label 'cpu_16'
    label 'mem_100M'
    label 'time_12'

    conda 'bioconda::gubbins=3.2.1'
    container 'quay.io/biocontainers/gubbins:3.2.1--py38pl5321h4c6a040_1'

    input:
    tuple path(msa), path(recombination_pred)

    output:
    path("masked_msa.aln"), emit: masked_msa

    script:
    """
    mask_gubbins_aln.py \\
        --missing-char N \\
        --aln ${msa} \\
        --gff ${recombination_pred} \\
        --out masked_msa.aln
    """
}
