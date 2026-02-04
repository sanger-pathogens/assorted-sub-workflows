process TRIMMOMATIC {
    tag "${meta.ID}"
    label 'mem_1'
    label 'time_1'
    label 'cpu_4'

    container "quay.io/biocontainers/trimmomatic:0.39--1"

    // publish only the gz version
    // need to be gzipped: need a process here to gzip the file.
    publishDir enabled: params.publish_trimmomatic_reads, mode: 'copy', pattern: "*_trimmed_1.fastq.gz", path: "${params.outdir}/${meta.ID}/preprocessing/trimmed_reads/"
    publishDir enabled: params.publish_trimmomatic_reads, mode: 'copy', pattern: "*_trimmed_2.fastq.gz", path: "${params.outdir}/${meta.ID}/preprocessing/trimmed_reads/"
    publishDir enabled: params.publish_trimmomatic_reads, mode: 'copy', pattern: "*_trimmed_unpaired_1.fastq.gz", path: "${params.outdir}/${meta.ID}/preprocessing/trimmed_reads/"
    publishDir enabled: params.publish_trimmomatic_reads, mode: 'copy', pattern: "*_trimmed_unpaired_2.fastq.gz", path: "${params.outdir}/${meta.ID}/preprocessing/trimmed_reads/"

    input:
    tuple val(meta), path(extracted_R1), path(extracted_R2)

    output:
    tuple val(meta), path(output_1), path(output_2), emit: trimmed_fastqs
    path(summary_file), emit: trimmomatic_stats

    script:
    output_1="${meta.ID}_trimmed_1.fastq"
    output_2="${meta.ID}_trimmed_2.fastq"
    output_1_unpaired="${meta.ID}_trimmed_unpaired_1.fastq"
    output_2_unpaired="${meta.ID}_trimmed_unpaired_2.fastq"
    summary_file="${meta.ID}_trimmomatic_summary.csv"

    gzip_cmd = params.publish_trimmomatic_reads ? """
    gzip -f ${output_1}
    gzip -f ${output_2}
    gzip -f ${output_1_unpaired}
    gzip -f ${output_2_unpaired}
    """ : ""

    """
    trimmomatic PE -phred33 -threads ${task.cpus} ${extracted_R1} ${extracted_R2} \
    ${output_1} ${output_1_unpaired} \
    ${output_2} ${output_2_unpaired} \
    ${params.trimmomatic_options} \
    -summary ${summary_file}

    ${gzip_cmd}
    """

}
