process QUAST {
    tag "${meta.ID}"
    label "cpu_1"
    label "mem_250M"
    label "time_30m"

    container  'quay.io/biocontainers/quast:5.0.2--py36pl5321hcac48a8_7'

    publishDir mode: 'copy', pattern: "${report_out}", path: "${params.outdir}/${meta.ID}/quast"

    input:
    tuple val(meta), path(consensus)

    output:
    tuple val(meta), path(report_out), emit: quast_out

    script:
    output = "${meta.ID}_assembly_stats"
    report_path = "${output}/transposed_report.tsv"
    report_txt = "${output}/report.txt"

    report_out = params.output_transposed ? report_path : report_txt

    """
    quast.py ${consensus} -o ${output} --no-html --no-plots
    """
}