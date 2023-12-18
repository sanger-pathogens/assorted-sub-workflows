process PICARD_MARKDUP {
    label 'cpu_2'
    label 'mem_2'
    label 'time_1'

    container 'quay.io/biocontainers/picard:3.1.1--hdfd78af_0'

    publishDir "${params.outdir}/${meta.id}/picard", enabled: params.keep_sorted_bam, mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(sorted_reads)

    output:
    tuple val(meta), path("${dedup_reads}"),  emit: dedup_reads

    script:
    dedup_reads = "${meta.id}_duplicates_removed.bam"
    """
    picard MarkDuplicates \
      -I $sorted_reads \
      -O $dedup_reads \
      -M marked_dup_metrics.txt \
      --REMOVE_DUPLICATES
    """
}

