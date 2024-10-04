process MULTIQC {
    label 'cpu_1'
    label 'mem_2'
    label 'time_30m'

    container 'quay.io/biocontainers/multiqc:1.19--pyhdfd78af_0'

    publishDir "${params.outdir}/multiqc/", pattern: "*.html", mode: 'copy', overwrite: true

    input:
    path('*') // we will likely need to collect + join into some mega channel unsure what to do here actually

    output:
    path("${date}-report.html"), emit: report
    path("*_data"), emit: data
    path("*_plots"), optional:true, emit: plots

    script:
    def out_report = "${outtag}-report-${date}.html"
    def custom_config = params.multiqc_config ? "--config ${params.multiqc_config}" : "" //add config if you supply one
    //workflow start is ugly 2024-02-29T12:01:26.233465Z split on T for time and take just the date
    date = "${workflow.start}".split('T')[0]
    """
    multiqc \
        -n ${out_report} \
        -f \
        ${custom_config} \
        .
    """
}