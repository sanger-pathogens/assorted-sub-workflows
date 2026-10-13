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
    label 'mem_8'
    label 'time_12'

    // minibwa 0.6 binary added to the samtools 1.22 biocontainer (params.minibwa_container).
    // The minibwa biocontainer has no samtools, so this used to write an uncompressed SAM
    // (median 12 GB, 3.6 TB over 300 study 8279 samples) for a separate SORT_TO_BAM process.
    container "${params.minibwa_container}"

    input:
    tuple val(meta), path(reads_1), path(reads_2), path(reference), path(bwa_index_files)

    output:
    tuple val(meta), path(mapped_reads_bam), emit: mapped_reads

    script:
    // -t threads, legacy bwa-mem-compatible CLI (`mem`, not the newer `map`) - same interface
    // as the original bwa mem call this replaces.
    // Nextflow runs this under bash -euo pipefail, and the explicit set keeps it that way if the
    // shell is overridden: without pipefail a minibwa crash would leave sort writing a truncated BAM.
    // Sort writes to .partial and is renamed only on success.
    mapped_reads_bam = "${meta.ID}_mapped.bam"
    sort_threads = Math.max(1, task.cpus.intdiv(2))
    """
    set -o pipefail
    minibwa mem -t ${task.cpus} ${reference} ${reads_1} ${reads_2} \\
      | samtools sort -@ ${sort_threads} -m 1G -T sort_tmp -o "${mapped_reads_bam}.partial" -
    mv "${mapped_reads_bam}.partial" "${mapped_reads_bam}"
    """
}
