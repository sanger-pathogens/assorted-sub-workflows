process BWA {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_12'
    label 'avx512'

    container 'quay.io/sangerpathogens/bwa:0.7.19'

    input:
    tuple val(meta), path(reads_1), path(reads_2), path(reference), path(bwa_index_files)

    output:
    tuple val(meta), path(mapped_reads),  emit: mapped_reads

    script:
    mapped_reads = "${meta.ID}_mapped.bam"
    // -v 1 for only errors, -M for picard compatibility
    """
    bwa mem -t ${task.cpus} ${reference} ${reads_1} ${reads_2} \
      | samtools view -@ ${task.cpus} -b - \
      | samtools sort -@ ${task.cpus} -o "${mapped_reads}"
    """

}

process BWA_INDEX {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_250M'
    label 'time_30m'
    label 'avx512'

    container 'quay.io/sangerpathogens/bwa:0.7.19'

    input:
    tuple val(meta), path(reference)

    output:
    tuple val(meta), path(reference), path("${reference}.*"),  emit: bwa_index

    script:
    """
    bwa index ${reference}
    """
}