process MODEL_FINDER {
    label 'cpu_2'
    label 'mem_16'
    label 'time_12'

    container 'quay.io/biocontainers/iqtree:2.3.6--hdbdd923_0'

    publishDir "${params.outdir}/tree/iqtree", mode: 'copy', overwrite: true, pattern: "${iqtree_log_liklihood_models}"

    input:
    path(msa)

    output:
    path(*)
    path(iqtree_log_liklihood_models), emit: inferred_models

    script:
    iqtree_log_liklihood_models = "${msa.baseName}.model.gz"
    """
    iqtree \\
        -s "${msa}" \\
        -m MF \\
        ${params.model_complexity_criterion}
    """
}