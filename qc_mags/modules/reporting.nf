process REPORT {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    publishDir mode: 'copy', path: "${params.outdir}/report/"

    container 'quay.io/sangerpathogens/pandas:2.2.1'

    input:
    tuple val(meta), path(pre_checkm2_report), path(pre_gunc_report), path(post_checkm2_report), path(post_gunc_report)

    output:
    tuple val(meta), path("*.csv"), emit: report

    script:
    command = "${projectDir}/assorted-sub-workflows/qc_mags/bin/report.py"
    """
    ${command} \\
        --pre_qc_checkm2 ${pre_checkm2_report}
        --pre_qc_gunc ${pre_gunc_report}
        --post_qc_checkm2 ${post_checkm2_report}
        --post_qc_gunc ${post_gunc_report}
        --output final_report.csv
    """
}