process MERGE_FASTQS {
    tag "${meta.id}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_queue_from_normal'
    
    input:
    tuple val(meta), path(first_read), path(second_read)

    output:
    tuple val(meta), path(merged_fastq), emit: merged_fastq

    script:
    merged_fastq="${meta.id}_merged.fastq.gz"
    """
    cat ${first_read} ${second_read} > ${meta.id}_merged.fastq.gz
    """
}
