process MULTIQC {
    label 'cpu_1'
    label 'mem_64'
    label 'time_30m'

    container 'quay.io/biocontainers/multiqc:1.19--pyhdfd78af_0'

    publishDir "${params.outdir}/multiqc/", mode: 'copy', overwrite: true

    input:
    path('*') // we will likely need to collect + join into some mega channel unsure what to do here actually

    output:
<<<<<<< HEAD
    path(output_report), emit: report
=======
    path(out_report), emit: report
>>>>>>> ccd4367 (Update multiqc module)
    path(output_data), emit: data
    path(output_plots), optional:true, emit: plots

    script:
    def custom_config = params.multiqc_config ? "--config ${params.multiqc_config}" : "" // add config if you supply one

<<<<<<< HEAD
    def date = "${workflow.start}".split('T')[0] // e.g. 2024-02-29T12:01:26.233465Z, so split on T to use only date
    output_report = "${date}-report.html"
=======
    date = "${workflow.start}".split('T')[0] // workflow start is ugly 2024-02-29T12:01:26.233465Z, so split on T to use only date
    out_report = "${date}-report.html"
>>>>>>> ccd4367 (Update multiqc module)
    output_data = "${date}_data.tar.gz"
    output_plots = "${date}_plots.tar.gz"

    """
    multiqc \\
<<<<<<< HEAD
        -n ${output_report} \\
=======
        -n ${out_report} \\
>>>>>>> ccd4367 (Update multiqc module)
        -f \\
        ${custom_config} \\
        .

    tar -czf ${output_data} \\
        -C *report_data \\
        --exclude 'multiqc_data.json' \\
        --transform='s,^./,${output_data.replaceFirst(/\.tar.gz$/, '')}/,' \\
        .

    if [[ -d ${date}-report_plots ]]; then
        tar -czf ${output_plots} \\
            -C *report_plots/svg  \\
            --transform='s,^./,${output_plots.replaceFirst(/\.tar.gz$/, '')}/,' \\
            .
    fi
    """
}