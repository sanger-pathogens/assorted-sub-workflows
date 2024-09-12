process KREPORT2MPA {
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    publishDir "${params.outdir}/${meta.ID}/bracken", mode: 'copy', overwrite: true, pattern: "*.mpa.txt"


    container 'quay.io/biocontainers/krakentools:1.2--pyh5e36f6f_0'

    input:
    tuple val(meta), path(kraken_style_bracken_report)

    output:
    tuple val(meta), path("${mpa_out}"),  emit: mpa_abundance_report

    script:
    mpa_out = "${meta.ID}_report_bracken_species.mpa.txt"
    """
    kreport2mpa.py -r ${kraken_style_bracken_report} -o ${mpa_out} --intermediate-ranks --display-header
    """
}

process GENERATE_ABUNDANCE_SUMMARY {
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    publishDir "${params.outdir}/abundance_summary", mode: 'copy', overwrite: true

    container 'quay.io/biocontainers/krakentools:1.2--pyh5e36f6f_0'

    input:
    path(mpa_abundance_reports)

    output:
    path("*summary_report.tsv"),  emit: abundance_summary

    script:
    summary_file = "bracken_summary_report.tsv"
    """
    # Combine in summary file
    combine_mpa.py -i *.mpa.txt -o "${summary_file}"
    # Fix header (get sample_names)
    sed -i 's;_kraken_sample_report_bracken_species.tsv;;g' "${summary_file}"
    """
}
