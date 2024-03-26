process COLLATE_FASTQ {
    label 'cpu_2'
    label 'mem_1'
    label "time_queue_from_${params.start_queue}"

    conda 'bioconda::samtools=1.17'
    container 'quay.io/biocontainers/samtools:1.17--hd87286a_2'

    publishDir path: { if ("${params.save_method}" == "nested") "${params.outdir}/${meta.ID}/${params.raw_reads_prefix}fastqs/" else "${params.outdir}/fastqs/" } , enabled: params.save_fastqs, mode: 'copy', overwrite: true, pattern: "*_1.fastq.gz", saveAs: { filename -> "${params.raw_reads_prefix}${forward_fastq}" }
    publishDir path: { if ("${params.save_method}" == "nested") "${params.outdir}/${meta.ID}/${params.raw_reads_prefix}fastqs/" else "${params.outdir}/fastqs/" } , enabled: params.save_fastqs, mode: 'copy', overwrite: true, pattern: "*_2.fastq.gz", saveAs: { filename -> "${params.raw_reads_prefix}${reverse_fastq}" }


    input:
    tuple val(meta), path(cram)

    output:
    tuple val(meta), path(forward_fastq), path(reverse_fastq), emit: fastq_channel
    val(meta), emit: metadata_channel
    path(cram), emit: remove_channel

    script:
    forward_fastq = "${meta.ID}_1.fastq.gz"
    reverse_fastq = "${meta.ID}_2.fastq.gz"
    """
    samtools collate -O \
    -f ${cram} \
    -@ ${task.cpus} \
    |
    samtools fastq -N \
        -1 ${forward_fastq} \
        -2 ${reverse_fastq} \
        -@ ${task.cpus}
    """
}

process COMBINE_FASTQ {
    label 'cpu_2'
    label 'mem_1'
    label 'time_12'

    conda 'bioconda::samtools=1.17'
    container 'quay.io/biocontainers/samtools:1.17--hd87286a_2'

    publishDir path: { if ("${params.save_method}" == "nested") "${params.outdir}/${meta.ID}/${params.raw_reads_prefix}fastqs/" else "${params.outdir}/fastqs/" } , enabled: params.save_fastqs, mode: 'copy', overwrite: true, pattern: "*_1.fastq.gz", saveAs: { filename -> "${params.raw_reads_prefix}${forward_fastq}" }
    publishDir path: { if ("${params.save_method}" == "nested") "${params.outdir}/${meta.ID}/${params.raw_reads_prefix}fastqs/" else "${params.outdir}/fastqs/" } , enabled: params.save_fastqs, mode: 'copy', overwrite: true, pattern: "*_2.fastq.gz", saveAs: { filename -> "${params.raw_reads_prefix}${reverse_fastq}" }
    
    input:
        tuple val(meta), val(read_1_paths), val(read_2_paths) // read_*_paths are strings containing multiple read file paths separated by a single space character

    output:
        tuple val(meta), path("${forward_fastq}.gz"), path("${reverse_fastq}.gz"), emit: fastq_channel

    script:
    forward_fastq = "${meta.ID}_1.fastq"
    reverse_fastq = "${meta.ID}_2.fastq"
    """
    zcat ${read_1} > ${forward_fastq}
    zcat ${read_2} > ${reverse_fastq}
    gzip ${forward_fastq}
    gzip ${reverse_fastq}
    """
}