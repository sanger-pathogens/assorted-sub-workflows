process REPORT {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    publishDir mode: 'copy', path: "${params.outdir}/report/"

    container 'quay.io/sangerpathogens/pandas:2.2.1'

    input:
    tuple(
        val(meta),
        path("pre_quast_report.tsv"),
        path("pre_checkm2_report.tsv"),
        path("pre_gunc_report.tsv"),
        path("gtdbtk_report.tsv"),
        path("post_quast_report.tsv"),
        path("post_checkm2_report.tsv"),
        path("post_gunc_report.tsv"),
        path(report_config)
    )

    output:
    tuple val(meta), path(final_report), emit: report

    script:
    command = "${projectDir}/assorted-sub-workflows/qc_mags/bin/report.py"
    final_report = "${meta.ID}_final_report.tsv"
    """
    ${command} \\
        --pre_qc_quast pre_quast_report.tsv \\
        --pre_qc_checkm2 pre_checkm2_report.tsv \\
        --pre_qc_gunc pre_gunc_report.tsv \\
        --gtdbtk gtdbtk_report.tsv \\
        --post_qc_quast post_quast_report.tsv \\
        --post_qc_checkm2 post_checkm2_report.tsv \\
        --post_qc_gunc post_gunc_report.tsv \\
        --config ${report_config} \\
        --output ${final_report}
    """
}

process SUMMARISE_CONTIG_FILTERING {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    publishDir mode: 'copy', path: "${params.outdir}/report/"

    container 'quay.io/sangerpathogens/pandas:2.2.1'

    input:
    tuple val(meta), path(final_report)

    output:
    tuple val(meta), path(filter_summary), emit: summary

    script:
    command = "${projectDir}/assorted-sub-workflows/qc_mags/bin/summarise_filtered_contigs.py"
    filter_summary = "${meta.ID}_filter_summary.tsv"

    """
    ${command} \\
        --report ${final_report} \\
        --min_contig ${params.min_contig} \\
        --output ${filter_summary}
    """
}