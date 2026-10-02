process COLLATE_FASTQ {
    label 'cpu_2'
    label "mem_${params.large_data ? 10 : 1}"
    label "time_queue_from_${params.large_data ? 'week' : 'normal'}"

    conda 'bioconda::samtools=1.17'
    container 'quay.io/biocontainers/samtools:1.17--hd87286a_2'

    publishDir path: { if ("${params.save_method}" == "nested") "${params.outdir}/${meta.ID}/${params.raw_reads_output_dir}/" else "${params.outdir}/${params.raw_reads_output_dir}/" } , enabled: params.save_fastqs, mode: 'copy', overwrite: true, pattern: "*_1.fastq.gz"
    publishDir path: { if ("${params.save_method}" == "nested") "${params.outdir}/${meta.ID}/${params.raw_reads_output_dir}/" else "${params.outdir}/${params.raw_reads_output_dir}/" } , enabled: params.save_fastqs, mode: 'copy', overwrite: true, pattern: "*_2.fastq.gz"


    input:
    val(meta)

    output:
    tuple val(meta), path(forward_fastq), path(reverse_fastq), emit: fastq_channel
    path(local_link), emit: files_to_remove

    script:
    input_file = file(meta.local_path)
    local_link = input_file.getName()
    forward_fastq = "${meta.ID}${params.raw_reads_suffix}_1.fastq.gz"
    reverse_fastq = "${meta.ID}${params.raw_reads_suffix}_2.fastq.gz"

    """
    ln -s ${input_file} ./${local_link}
    samtools collate -O \
    -f ${local_link} \
    -@ ${task.cpus} \
    |
    samtools fastq -N \
        -1 ${forward_fastq} \
        -2 ${reverse_fastq} \
        -@ ${task.cpus}
    """
}
