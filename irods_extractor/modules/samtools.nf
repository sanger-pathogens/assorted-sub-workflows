process COLLATE_FASTQ {
    label 'cpu_2'
    label "mem_${params.large_data ? 10 : 1}"
    label "time_queue_from_${params.large_data ? 'week' : 'normal'}"

    conda 'bioconda::samtools=1.17'
    container 'quay.io/biocontainers/samtools:1.17--hd87286a_2'

    publishDir path: { if ("${params.save_method}" == "nested") "${params.outdir}/${meta.ID}/${params.raw_reads_prefix}fastqs/" else "${params.outdir}/fastqs/" } , enabled: params.save_fastqs, mode: 'copy', overwrite: true, pattern: "*_1.fastq.gz", saveAs: { filename -> "${params.raw_reads_prefix}${forward_fastq}" }
    publishDir path: { if ("${params.save_method}" == "nested") "${params.outdir}/${meta.ID}/${params.raw_reads_prefix}fastqs/" else "${params.outdir}/fastqs/" } , enabled: params.save_fastqs, mode: 'copy', overwrite: true, pattern: "*_2.fastq.gz", saveAs: { filename -> "${params.raw_reads_prefix}${reverse_fastq}" }


    input:
    val(meta)

    output:
    tuple val(meta), path(forward_fastq), path(reverse_fastq), emit: fastq_channel

    script:
    forward_fastq = "${meta.ID}_1.fastq.gz"
    reverse_fastq = "${meta.ID}_2.fastq.gz"
    """
    samtools collate -O \
    -f ${meta.local_path} \
    -@ ${task.cpus} \
    |
    samtools fastq -N \
        -1 ${forward_fastq} \
        -2 ${reverse_fastq} \
        -@ ${task.cpus}
    """
}
