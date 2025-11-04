process QUAST {
    tag "${meta.ID}"
    label "cpu_1"
    label "mem_1"
    label "time_30m"

    container   'quay.io/biocontainers/quast:5.3.0--py39pl5321heaaa4ec_0'

    publishDir mode: 'copy', path: "${params.outdir}/quast/"

    input:
    tuple val(meta), path(fastas, stageAs: "fastas/*")

    output:
    tuple val(meta), path(quast_report), emit: results

    script:
    quast_report = "${meta.ID}_quast_report.tsv"
    """
    quast.py fastas/*.${params.fasta_ext} -o quast --no-html --no-plots --min-contig ${params.min_contig}

    mv quast/transposed_report.tsv ${quast_report}
    """
}

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
    command = "${projectDir}/assorted-sub-workflows/qc_isolates/bin/quast_summary.py"
    report_tsv = "${meta.ID}_quast_summary.tsv"
    """
    ${command} \\
        --input ${quast_report} \\
        --output ${report_tsv}
    """
}