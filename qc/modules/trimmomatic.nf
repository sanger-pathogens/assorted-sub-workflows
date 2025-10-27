// REQUIRED
params.trimmomatic_options = "ILLUMINACLIP:${params.adapter_fasta}:2:10:7:1 CROP:151 SLIDINGWINDOW:4:20 MINLEN:70"
params.trimmomatic_threads = 4

process TRIMMOMATIC {
    tag "${meta.id}"
    label 'mem_1'
    label 'time_1'
    cpus params.trimmomatic_threads

    container "quay.io/biocontainers/trimmomatic:0.39--1"

    // publish only the gz version
    publishDir enabled: params.debug_preproc_output, mode: 'copy', failOnError: true, pattern: "${output_1_gz}", path: "${params.results_dir}/${meta.id}/preprocessing/trimmed_reads/"
    publishDir enabled: params.debug_preproc_output, mode: 'copy', failOnError: true, pattern: "${output_2_gz}", path: "${params.results_dir}/${meta.id}/preprocessing/trimmed_reads/"

    input:
    tuple val(meta), path(extracted_R1), path(extracted_R2)

    output:
    tuple val(meta), path(output_1), path(output_2), emit: paired_channel
    tuple val(meta), path(output_1_gz), path(output_2_gz)

    script:
    output_1="${meta.id}_trimmed_1.fastq"
    output_2="${meta.id}_trimmed_2.fastq"
    output_1_gz = "${output_1}.gz"
    output_2_gz = "${output_2}.gz"
    output_1_unpaired="${meta.id}_trimmed_unpaired_1.fastq"
    output_2_unpaired="${meta.id}_trimmed_unpaired_2.fastq"
    """
    trimmomatic PE -phred33 -threads ${params.trimmomatic_threads} ${extracted_R1} ${extracted_R2} \
    ${output_1} ${output_1_unpaired} \
    ${output_2} ${output_2_unpaired} \
    ${params.trimmomatic_options}
    gzip -c ${output_1} > ${output_1}.tmp.gz
    gzip -c ${output_2} > ${output_2}.tmp.gz
    mv ${output_1}.tmp.gz ${output_1_gz}
    mv ${output_2}.tmp.gz ${output_2_gz}
    """
}
