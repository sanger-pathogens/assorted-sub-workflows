process GENERATE_ABUNDANCE_SUMMARY {
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    publishDir "${params.outdir}/abundance_summary", mode: 'copy', overwrite: true

    container 'quay.io/biocontainers/krakentools:1.2--pyh5e36f6f_0'

    input:
    path(kraken_style_bracken_reports)

    output:
    path("*summary_report.tsv"),  emit: abundance_summary

    script:
    summary_file = "summary_report.tsv"
    """
    mpa_reports="mpa_reports"
    # Convert bracken reports to mpa format
    mkdir -p "\${mpa_reports}"
    for report in *.tsv; do
        kreport2mpa.py -r "\${report}" -o "\${mpa_reports}/\${report}" --intermediate-ranks --display-header
    done
    # Combine in summary file
    combine_mpa.py -i \${mpa_reports}/* -o "${summary_file}"
    # Fix header (get sample_names)
    sed -i 's;_kraken_sample_report_bracken_[^.]*\\.tsv;;g' "${summary_file}"
    """
}
