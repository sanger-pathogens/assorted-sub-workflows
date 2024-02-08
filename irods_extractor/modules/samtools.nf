process COLLATE_CRAM {
    label 'cpu_2'
    label 'mem_1'
    label 'time_12'

    conda 'bioconda::samtools=1.17'
    container 'quay.io/biocontainers/samtools:1.17--hd87286a_2'
    input:
    tuple val(meta), path(cram)

    output:
    tuple val(meta), path(bam), path(cram), emit: bam_channel

    script:
    bam = "${meta.ID}.bam"
    """
    samtools collate -u -o ${bam} -f ${cram} -@ ${task.cpus}
    """
}

process FASTQ_FROM_COLLATED_BAM {
    label 'cpu_2'
    label 'mem_1'
    label 'time_12'

    conda 'bioconda::samtools=1.17'
    container 'quay.io/biocontainers/samtools:1.17--hd87286a_2'

    publishDir path: { if ("${params.save_method}" == "nested") "${params.outdir}/${meta.ID}/${params.raw_reads_prefix}fastqs/" else "${params.outdir}/fastqs/" } , enabled: params.save_fastqs, mode: 'copy', overwrite: true, pattern: "*_1.fastq.gz", saveAs: { filename -> "${params.raw_reads_prefix}${forward_fastq}" }
    publishDir path: { if ("${params.save_method}" == "nested") "${params.outdir}/${meta.ID}/${params.raw_reads_prefix}fastqs/" else "${params.outdir}/fastqs/" } , enabled: params.save_fastqs, mode: 'copy', overwrite: true, pattern: "*_2.fastq.gz", saveAs: { filename -> "${params.raw_reads_prefix}${reverse_fastq}" }


    input:
    tuple val(meta), path(bam), path(cram)

    output:
    tuple val(meta), path(forward_fastq), path(reverse_fastq), emit: fastq_channel
    val(meta), emit: metadata_channel
    tuple path(bam), path(cram), emit: remove_channel

    script:
    forward_fastq = "${meta.ID}_1.fastq.gz"
    reverse_fastq = "${meta.ID}_2.fastq.gz"
    """
    samtools fastq -N \
        -1 ${forward_fastq} \
        -2 ${reverse_fastq} \
        -@ ${task.cpus} \
        ${bam}
    """
}