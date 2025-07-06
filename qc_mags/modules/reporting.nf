process REPORT {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    publishDir mode: 'copy', path: "${params.outdir}/report/"

    container 'quay.io/sangerpathogens/python-curl:3.11'

    input:
    tuple val(meta), path(merged_csv), path(contigs)

    output:
    tuple val(meta), path("*.csv"), emit: report

    script:
    command = "${projectDir}/bin/report.py"
    """
    ${command} ${merged_csv}
    """
}