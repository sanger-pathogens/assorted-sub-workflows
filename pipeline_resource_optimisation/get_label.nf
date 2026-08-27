process GET_LABEL {
    label 'cpu_1'
    label 'mem_500M'
    label 'time_12'

    container "quay.io/sangerpathogens/pipeline_resource_optimisation:test_6"

    output:
    path("label.txt"), emit: label_file

    script:
    """
    modelling \
        --mode get_label \
        -n METASPADES_GradientBoosting \
        -tag 0.0.0 \
        --predict \
        -d /lustre/scratch124/pam/teams/team230/tm22/tickets/PAT-3483/sra_test/DRR671405.json \
        --adjust-value mae \
        --adjust-scale 1.5 \
    > label.txt
    """
}



