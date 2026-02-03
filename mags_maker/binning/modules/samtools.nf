process INDEX {
    label 'cpu_1'
    label 'mem_100M'
    label 'time_30m'

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
