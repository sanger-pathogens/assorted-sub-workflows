process GET_LABEL {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_500M'
    label 'time_12'

    container "quay.io/sangerpathogens/pipeline_resource_optimisation:test_7"

    input:
    tuple val(meta), val(features_json)

    output:
    tuple val(meta), stdout, emit: label

    script:
    // features_json is a single JSON object holding one row of the model x table.
    // The predicted value is left on stdout so it can go straight into the meta map -
    // anything else the tool has to say must go to stderr to keep the value clean.
    """
    cat <<'END_FEATURES_JSON' > features.json
    ${features_json}
    END_FEATURES_JSON

    modelling \\
        --mode get_label \\
        -n ${params.model} \\
        -tag ${params.model_tag} \\
        --predict \\
        -d features.json \\
        --adjust-value ${params.model_adjust_value} \\
        --adjust-scale ${params.model_adjust_scale}
    """
}
