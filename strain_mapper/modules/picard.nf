process PICARD_MARKDUP {
    label 'cpu_2'
    label 'mem_2'
    label 'time_1'

    conda 'bioconda::picard=3.1.1'
    container 'quay.io/biocontainers/picard:3.1.1--hdfd78af_0'

    publishDir "${params.outdir}/${meta.ID}/picard", enabled: params.keep_sorted_bam, mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(sorted_reads)

    output:
    tuple val(meta), path("${dedup_reads}"),  emit: dedup_reads

    script:
    dedup_reads = "${meta.ID}_duplicates_removed.bam"
    """
    picard MarkDuplicates \
      -I $sorted_reads \
      -O $dedup_reads \
      -M marked_dup_metrics.txt \
      --REMOVE_DUPLICATES
    """
}

process PICARD_ADD_READGROUP {
    label 'cpu_2'
    label 'mem_2'
    label 'time_1'

    container 'quay.io/biocontainers/picard:3.1.1--hdfd78af_0'

    publishDir "${params.outdir}/${meta.id}/picard", enabled: params.keep_sorted_bam, mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(sorted_reads)

    output:
    tuple val(meta), path("${rg_added_reads}"),  emit: rg_added_reads

    script:
    rg_added_reads = "${sorted_reads.baseName}_rg_added.bam"
    """
    picard AddOrReplaceReadGroups \
      I=${sorted_reads} \
      O=${rg_added_reads} \
      RGLB=NA \
      RGPL=NA \
      RGPU=NA \
      RGSM=${meta.id}
    """
}
