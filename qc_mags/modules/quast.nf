process QUAST {
    tag "${meta.ID}"
    label "cpu_1"
    label "mem_8"
    label "time_30m"

    container   'quay.io/biocontainers/quast:5.3.0--py39pl5321heaaa4ec_0'

    publishDir mode: 'copy', path: "${params.outdir}/quast/"

    input:
    tuple val(meta), path(fastas, stageAs: "fastas/*")

    output:
    tuple val(meta), path(report_tsv), emit: results

    script:
    report_tsv = "${meta.ID}_quast_summary.tsv"
    """
    quast.py fastas/*.${params.fasta_ext} -o quast --no-html --no-plots

    mv quast/transposed_report.tsv ${report_tsv}
    """
}