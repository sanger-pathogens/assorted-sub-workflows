def getQuastThresholds() {
    // Uses QUAST default contig thresholds and min_contig if value not already in the defaults
    // Threshold (x) used for reporting # contigs > x and total length of contigs > x
    def quast_contig_thresholds = [0,1000,5000,10000,25000,50000]
    if (params.min_contig) {
        def min_contig = params.min_contig as Integer
        if (!(min_contig in quast_contig_thresholds)) {
            quast_contig_thresholds += min_contig
            quast_contig_thresholds = quast_contig_thresholds.sort()
        }
    }
    // Outputs as a comma-separated list (string) ready as input to quast command param contig-thresholds
    def thresholds_str = quast_contig_thresholds.join(',')
    return thresholds_str
}

process QUAST {
    tag "${meta.ID}"
    label "cpu_1"
    label "mem_1"
    label "time_30m"

    container   'quay.io/biocontainers/quast:5.3.0--py39pl5321heaaa4ec_0'

    publishDir mode: 'copy', path: "${params.outdir}/${qc_stage}/quast/", enabled: !(params.skip_raw_reports)

    input:
    tuple val(meta), path(fastas, stageAs: "fastas/*")
    val(qc_stage)
    val(thresholds_str)

    output:
    tuple val(meta), path(quast_report), emit: results

    script:
    quast_report = "${meta.ID}_quast_report.tsv"
    """
    quast.py fastas/*.${params.fasta_ext} -o quast --no-html --no-plots --contig-thresholds ${thresholds_str} --threads ${task.cpus}

    mv quast/transposed_report.tsv ${quast_report}
    """
}

process QUAST_SUMMARY {
    tag "${meta.ID}"
    label "cpu_1"
    label "mem_1"
    label "time_1"

    container 'quay.io/sangerpathogens/pandas:2.2.1'

    publishDir mode: 'copy', path: "${params.outdir}/${qc_stage}/quast_summary/", enabled: !(params.skip_raw_reports)

    input:
    tuple val(meta), path(quast_report)
    val(qc_stage)

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