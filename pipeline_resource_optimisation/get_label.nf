process GET_LABEL {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_250M'
    label 'time_12'

    input:
    tuple val(meta), path(model), path(predictors)

    output:
    tuple val(meta), env(RESOURCE_LABEL)

    script:
    def command = "${projectDir}/assorted-sub-workflows/pipeline_resource_optimisation/get_label.py"
    """
    RESOURCE_LABEL=\$(
        ${command} \
            --model ${model} \
            --predictors ${predictors}
    )
    """
}