process PICARD_MARKDUP {
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'

    container 'quay.io/biocontainers/picard:3.1.1--hdfd78af_0'

    input:
    tuple val(meta), path("${sorted_reads}")

    output:
    tuple val(meta), path("${dedup_reads}"),  emit: dedup_reads

    script:
    mpileup_file = "${meta.id}_marked_duplicates.bam"
    """
    picard MarkDuplicates \
      --REMOVE_DUPLICATES \
      I=$sorted_reads \
      O=$mpileup_file \
      M=marked_dup_metrics.txt
    """
}

