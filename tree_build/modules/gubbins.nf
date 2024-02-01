process GUBBINS{
    label 'cpu_16'
    label 'mem_1'
    label 'time_12'

    conda 'bioconda::gubbins=3.2.1'
    container 'quay.io/biocontainers/gubbins:3.2.1--py38pl5321h4c6a040_1'

    publishDir "${params.outdir}/gubbins", mode: 'copy', overwrite: true, pattern: "${gubprefix}.*"

    input:
    path(msa)

    output:
    tuple path(msa), path(recpredgff), emit: recpredgff_ch

    script:
    gubprefix="gubbins_out"
    recpredgff="${gubprefix}.recombination_predictions.gff"
    recfreetree="${gubprefix}.node_labelled.final_tree.tre"
    """
    run_gubbins.py --prefix ${gubprefix} \
    --filter-percentage ${params.gubbins_filter_percentage} \
    --threads ${params.gubbins_threads} \
    ${msa}
    """
}

process GUBBINS_MASK{
    label 'cpu_16'
    label 'mem_100M'
    label 'time_12'

    conda 'bioconda::gubbins=3.2.1'
    // need to check if equivalent to below: container '/software/pathogen/images/gubbins-3.2.1.simg'
    container 'quay.io/biocontainers/gubbins:3.2.1--py38pl5321h4c6a040_1'


    input:
    tuple path(msa), path(recpredgff)

    output:
    path("masked_msa.aln"), emit: masked_msa

    script:
    """
    mask_gubbins_aln.py --missing-char N --aln ${msa} --gff ${recpredgff} --out masked_msa.aln
    """
}
