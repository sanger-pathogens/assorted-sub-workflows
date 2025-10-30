process COMPRESS_READS {
    tag "${meta.ID}"
    label 'mem_4'
    label 'cpu_1'
    label 'time_queue_from_small_slow2'

    publishDir "${params.outdir}/${meta.ID}/preprocessing/", mode: "copy", enabled: params.publish_clean_reads

    input:
    tuple val(meta), path(read_1), path(read_2)

    output:
    tuple val(meta), path("${meta.ID}_pre_processed_1.fq.gz"), path("${meta.ID}_pre_processed_2.fq.gz")

    script:
    """
    gzip -c ${read_1} > ${read_1}.tmp.gz
    gzip -c ${read_2} > ${read_2}.tmp.gz
    mv ${read_1}.tmp.gz ${meta.ID}_pre_processed_1.fq.gz
    mv ${read_2}.tmp.gz ${meta.ID}_pre_processed_2.fq.gz
    """
}
