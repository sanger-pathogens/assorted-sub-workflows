process TRIMMOMATIC {
    tag "${meta.ID}"
    label 'mem_1'
    label 'time_1'
    label 'cpu_4'
    // cpus params.trimmomatic_threads

    container "quay.io/biocontainers/trimmomatic:0.39--1"

    // publish only the gz version
    publishDir enabled: params.debug_preproc_output, mode: 'copy', failOnError: true, pattern: "${output_1_gz}", path: "${params.outdir}/${meta.ID}/preprocessing/trimmed_reads/"
    publishDir enabled: params.debug_preproc_output, mode: 'copy', failOnError: true, pattern: "${output_2_gz}", path: "${params.outdir}/${meta.ID}/preprocessing/trimmed_reads/"

    input:
    tuple val(meta), path(extracted_R1), path(extracted_R2)

    output:
    tuple val(meta), path(output_1), path(output_2), emit: trimmed_fastqs

    script:
    output_1="${meta.ID}_trimmed_1.fastq"
    output_2="${meta.ID}_trimmed_2.fastq"
    output_1_unpaired="${meta.ID}_trimmed_unpaired_1.fastq"
    output_2_unpaired="${meta.ID}_trimmed_unpaired_2.fastq"
    """
    trimmomatic PE -phred33 -threads ${task.cpus} ${extracted_R1} ${extracted_R2} \
    ${output_1} ${output_1_unpaired} \
    ${output_2} ${output_2_unpaired} \
    ${params.trimmomatic_options}
    """
}
