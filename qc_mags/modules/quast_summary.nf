process QUAST_SUMMARY {
    tag "${meta.ID}"
    label "cpu_1"
    label "mem_1"
    label "time_1"

    container 'quay.io/sangerpathogens/pandas:2.2.1'

    publishDir mode: 'copy', path: "${params.outdir}/quast_summary/"

    input:
    tuple val(meta), path(quast_report)

    output:
    tuple val(meta), path(report_tsv), emit: results

    script:
    command = "${projectDir}/assorted-sub-workflows/qc_mags/bin/quast_summary.py"
    report_tsv = "${meta.ID}_quast_summary.tsv"
    """
    ${command} \\
        --input ${quast_report} \\
        --output ${report_tsv}
    """
}