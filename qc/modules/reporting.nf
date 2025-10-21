process REPORT {
    label 'cpu_1'
    label 'mem_2'
    label 'time_30m'

    publishDir "${params.outdir}/qc_pass_fail_summary/", pattern: "sample_pass_fail_qc_summary.tsv", mode: 'copy', overwrite: true

    input:
    val(pass_fail_channel)

    output:
    path("sample_pass_fail_qc_summary.tsv")

    script:

    def header = ['Sample_ID', 'FastQC_Status']
        if (params.bracken_profile) header << 'Kraken_Braken_Profile_Status'
        if (params.sylph_profile) header << 'Sylph_Profile_Status'
    header = header.join('\t')

    def rows = pass_fail_channel.collect { row ->
        row.join('\t')
    }.join('\n')

    """
    echo -e "${header}\n${rows}" > sample_pass_fail_qc_summary.tsv
    """
}