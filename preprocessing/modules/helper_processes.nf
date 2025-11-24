process COMPRESS_READS {
    tag "${meta.ID}"
    label 'mem_4'
    label 'cpu_1'
    label 'time_queue_from_small_slow2'

    publishDir "${params.outdir}/${meta.ID}/preprocessing/", mode: "copy", enabled: params.publish_clean_reads

    input:
    tuple val(meta), path(read_1), path(read_2)

    output:
    tuple val(meta), path("${meta.ID}_preprocessed_1.fastq.gz"), path("${meta.ID}_preprocessed_2.fastq.gz")

    script:
    """
    gzip -c ${read_1} > ${read_1}.tmp.gz
    gzip -c ${read_2} > ${read_2}.tmp.gz
    mv ${read_1}.tmp.gz ${meta.ID}_preprocessed_1.fastq.gz
    mv ${read_2}.tmp.gz ${meta.ID}_preprocessed_2.fastq.gz
    """
}

process DECOMPRESS_READS {
    tag "${meta.ID}"
    label 'mem_4'
    label 'cpu_1'
    label 'time_queue_from_small_slow2'

    input:
    tuple val(meta), path(read_1), path(read_2)

    output:
    tuple val(meta), path("${meta.ID}_1.fastq"), path("${meta.ID}_2.fastq")

    script:
    """
    gunzip -c ${read_1} > ${meta.ID}_1.fastq
    gunzip -c ${read_2} > ${meta.ID}_2.fastq
    """
}
