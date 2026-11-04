process FASTP {
    tag "${meta.ID}"
    label 'cpu_4'
    label 'mem_4'
    label 'time_queue_from_normal'

    container 'quay.io/biocontainers/fastp:0.23.4--hadf994f_2'

    publishDir path: { if ("${params.save_method}" == "nested") "${params.outdir}/${meta.ID}/${params.preprocessing_output_dir}/" else "${params.outdir}/${params.preprocessing_output_dir}/" }, pattern: "*.fastq.gz", mode: 'copy', overwrite: true, enabled: params.publish_clean_reads
    publishDir "${params.outdir}/${meta.ID}/fastp/", pattern: "*.fastp.{json,html}", mode: 'copy', overwrite: true, enabled: params.save_individual_fastp_reports

    input:
    tuple val(meta), path(read_1), path(read_2)

    output:
    tuple val(meta), path("${preproc_read_stem}_1.fastq.gz"), path("${preproc_read_stem}_2.fastq.gz"), emit: reads
    path("${meta.ID}.fastp.json"), emit: json
    path("${meta.ID}.fastp.html"), emit: html

    script:
    preproc_read_stem = "${meta.ID}${params.preproc_read_suffix}"

    def poly_g = params.fastp_trim_poly_g ? "-g" : "-G"
    """
    fastp -w ${task.cpus} -i ${read_1} -I ${read_2} \\
        -o ${preproc_read_stem}_1.fastq.gz -O ${preproc_read_stem}_2.fastq.gz \\
        ${params.fastp_options} ${poly_g} \\
        -j ${meta.ID}.fastp.json -h ${meta.ID}.fastp.html
    """
}

process COLLATE_FASTP_STATS {
    tag "all samples"
    label 'cpu_1'
    label 'mem_1'
    label 'time_30m'

    container 'quay.io/sangerpathogens/python-curl:3.11'

    publishDir "${params.outdir}/preprocessing_summary_stats", mode: 'copy', overwrite: true, pattern: "*_statistics.tsv"

    input:
    path(fastp_jsons)

    output:
    path("fastp_statistics.tsv"), emit: fastp_stats

    script:
    command = "${projectDir}/assorted-sub-workflows/qc/bin/fastp_stats.py"
    """
    ${command} ${fastp_jsons} > fastp_statistics.tsv
    """
}
