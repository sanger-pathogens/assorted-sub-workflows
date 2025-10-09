process REPORT {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    publishDir mode: 'copy', path: "${params.outdir}/report/"

    container 'quay.io/sangerpathogens/pandas:2.2.1'

    input:
    tuple val(meta), path("checkm2_report.tsv"), path("gunc_report.tsv"), path("gtdbtk_report.tsv"), path("quast_summary_report.tsv"), path(report_config)

    output:
    tuple val(meta), path(final_report), emit: report

    script:
    command = "${projectDir}/assorted-sub-workflows/qc_isolates/bin/report.py"
    final_report = "${meta.ID}_final_report.tsv"
    """
    ${command} \\
        --checkm2 checkm2_report.tsv \\
        --gunc gunc_report.tsv \\
        --gtdbtk gtdbtk_report.tsv \\
        --quast quast_summary_report.tsv \\
        --config ${report_config} \\
        --output ${final_report}
    """
}