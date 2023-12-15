process PICARD_MARKDUP {
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'

    container 'quay.io/biocontainers/picard:3.1.1--hdfd78af_0'

    if (params.keep_sorted_bam){ publishDir "${params.outdir}/${meta.id}/picard", mode: 'copy', overwrite: true }

    input:
    tuple val(meta), path(sorted_reads)

    output:
    tuple val(meta), path("${dedup_reads}"),  emit: dedup_reads

    script:
    dedup_bam = "${meta.id}_duplicates_removed.bam"
    """
    picard MarkDuplicates \
      -I $sorted_reads \
      -O $dedup_bam \
      -M marked_dup_metrics.txt \
      --REMOVE_DUPLICATES
    """
}

