process MINIBWA_INDEX {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_2'
    label 'time_12'

    container 'quay.io/sangerpathogens/minibwa-samtools:0.7-1.23'

    input:
    tuple val(meta), path(reference)

    output:
    tuple val(meta), path(reference), path("${reference}.*"),  emit: bwa_index

    script:
    """
    minibwa index ${reference}
    """
}

process MINIBWA {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_2'
    label 'time_12'

    container 'quay.io/sangerpathogens/minibwa-samtools:0.7-1.23'

    input:
    tuple val(meta), path(reads_1), path(reads_2), path(reference), path(bwa_index_files)

    output:
    tuple val(meta), path(mapped_bam), emit: mapped_bam

    script:
    mapped_bam = "${meta.ID}_mapped.bam"
    """
    minibwa mem -t ${task.cpus} ${reference} ${reads_1} ${reads_2} \\
      | samtools view -@ ${task.cpus} -b - \\
      | samtools sort -@ ${task.cpus} -o "${mapped_bam}"
    """
}
