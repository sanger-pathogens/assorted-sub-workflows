process PICARD_MARKDUP {
    label 'cpu_2'
    label 'mem_2'
    label 'time_1'

    conda 'bioconda::picard=3.1.1'
    container 'quay.io/biocontainers/picard:3.1.1--hdfd78af_0'

    publishDir "${params.outdir}/${meta.ID}/picard", enabled: params.keep_dedup_bam, mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(sorted_reads_bam), path(sorted_reads_bai)

    output:
    tuple val(meta), path("${dedup_reads_bam}"),  emit: dedup_reads

    script:
    dedup_reads_bam = "${meta.ID}_duplicates_removed.bam"
    """
    picard MarkDuplicates \
      -I ${sorted_reads_bam} \
      -O ${dedup_reads_bam} \
      -M marked_dup_metrics.txt \
      --REMOVE_DUPLICATES
    """
}

