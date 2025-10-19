process REPORT {
    label 'cpu_1'
    label 'mem_2'
    label 'time_30m'

    publishDir "${params.outdir}/qc_pass_fail_summary/", pattern: "sample_pass_fail_qc_summary.tsv", mode: 'copy', overwrite: true, enabled: true

    input:
    val results

    output:
    path("sample_pass_fail_qc_summary.tsv")

    script:
    def header = "Sample_ID\tFastQC_Status\tKraken_Braken_Profile_Status\tSylph_Profile_Status"
    def rows = results.collect { row ->
        def (sample_meta, fastqc, kraken, sylph) = row
        "${sample_meta.ID}\t${fastqc}\t${kraken}\t${sylph}"
    }.join("\n")

    """
    echo -e "${header}\n${rows}" > sample_pass_fail_qc_summary.tsv
    """
}