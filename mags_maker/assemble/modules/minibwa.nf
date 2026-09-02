process MINIBWA_INDEX {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_250M'
    label 'time_12'

    container 'quay.io/biocontainers/minibwa:0.6--hab16a5f_0'

    input:
    tuple val(meta), path(reference)

    output:
    tuple val(meta), path(reference), path("${reference}.*"),  emit: bwa_index

    script:
    """
    minibwa index ${reference}
    """
}

process MINIBWA {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_12'

    container 'quay.io/biocontainers/minibwa:0.6--hab16a5f_0'

    input:
    tuple val(meta), path(reads_1), path(reads_2), path(reference), path(bwa_index_files)

    output:
    tuple val(meta), path(mapped_sam), emit: mapped_sam

    script:
    // minibwa's container ships only minibwa itself (no samtools), so SAM->sorted BAM
    // conversion is a separate process (SORT_TO_BAM in samtools.nf) - see MINIBWA_WF below.
    // -t threads, legacy bwa-mem-compatible CLI (`mem`, not the newer `map`) - same interface
    // as the original bwa mem call this replaces.
    mapped_sam = "${meta.ID}_mapped.sam"
    """
    minibwa mem -t ${task.cpus} ${reference} ${reads_1} ${reads_2} > ${mapped_sam}
    """
}
