process MULTIQC {
    label 'cpu_1'
    label 'mem_64'
    label 'time_30m'

    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'

    publishDir "${params.outdir}/multiqc/", mode: 'copy', overwrite: true

    input:
    path('*') // we will likely need to collect + join into some mega channel unsure what to do here actually
    path(multiqc_config)
    val(report_name_prefix)

    output:
    path(out_report), emit: report
    path(output_data), emit: data
    path(output_plots), optional:true, emit: plots

    script:
    def custom_config = multiqc_config ? "--config ${multiqc_config}" : "" // add config if you supply one

    date = "${workflow.start}".split('T')[0] // workflow start is ugly 2024-02-29T12:01:26.233465Z, so split on T to use only date
    out_report_name_base = "${date}-${report_name_prefix}"
    out_report = "${out_report_name_base}.html"
    output_data_base = "${out_report_name_base}_data"
    output_data = "${output_data_base}.tar.gz"
    output_plots_base = "${out_report_name_base}_plots"
    output_plots = "${output_plots_base}.tar.gz"

    """
    multiqc \\
        -n ${out_report} \\
        -f \\
        ${custom_config} \\
        .

    tar -czf ${output_data} \\
        -C ${output_data_base} \\
        --exclude 'multiqc_data.json' \\
        --transform='s,^./,${output_data_base}/,' \\
        .

    if [[ -d ${output_plots_base} ]]; then
        tar -czf ${output_plots} \\
            -C ${output_plots_base}/svg  \\
            --transform='s,^./,${output_plots_base}/,' \\
            .
    fi
    """
}