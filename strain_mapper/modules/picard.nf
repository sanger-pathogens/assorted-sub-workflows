process PICARD_MARKDUP {
    label 'cpu_2'
    label 'mem_2'
    label 'time_1'

    conda 'bioconda::picard=3.1.1'
    // TO DO when module installed on farm: 
    // container "${ singularity.enabled ? '/software/pathogen/images/picard:3.1.1--hdfd78af_0.simg' : 'quay.io/biocontainers/picard:3.1.1--hdfd78af_0' }"
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

