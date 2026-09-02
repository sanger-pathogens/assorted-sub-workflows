process INDEX {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_100M'
    label 'time_12'

    container 'quay.io/biocontainers/samtools:1.22--h96c455f_0'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path(bam), path("*.bam.*"), emit: bam_plus_index

    script:
    mapped_reads_bam = "${meta.ID}.bam"
    """
    samtools index ${bam}
    """
}

process SORT_TO_BAM {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_12'

    container 'quay.io/biocontainers/samtools:1.22--h96c455f_0'

    input:
    tuple val(meta), path(sam)

    output:
    tuple val(meta), path(mapped_reads_bam), emit: mapped_reads

    script:
    // minibwa's container ships only minibwa itself, no samtools - this is the SAM->sorted
    // BAM conversion step that used to be piped inline within bwa.nf's BWA process.
    mapped_reads_bam = "${meta.ID}_mapped.bam"
    """
    samtools view -@ ${task.cpus} -b ${sam} \\
      | samtools sort -@ ${task.cpus} -o "${mapped_reads_bam}"
    """
}
