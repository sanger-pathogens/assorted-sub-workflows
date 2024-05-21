process BWA {
    label 'cpu_1'
    label 'mem_8'
    label 'time_12'

    conda "bioconda::bwa=0.7.17"
    container 'quay.io/biocontainers/bwa:0.7.17--he4a0461_11'

    input:
    tuple val(meta), path(reads_1), path(reads_2)
    tuple path(reference), path(bwa_index_files)

    output:
    tuple val(meta), path("${meta.ID}.sam"),  emit: mapped_reads

    script:
    mapped_reads = "${meta.ID}.sam"
    // -v 1 for only errors -M for picard compatibility -a output all alignements
    """
    bwa mem -v 1 -M -a -t ${task.cpus}  ${reference} ${reads_1} ${reads_2} > ${meta.ID}.sam
    """
}

process BWA_INDEX {
    label 'cpu_4'
    label 'mem_8'
    label 'time_12'

    publishDir "${params.outdir}/bwa", mode: 'copy', overwrite: true

    conda "bioconda::bwa=0.7.17"
    container 'quay.io/biocontainers/bwa:0.7.17--he4a0461_11'

    input:
    path(reference)

    output:
    tuple path(reference), path("${reference}.*"),  emit: bwa_index

    script:
    """
    bwa index  ${reference}
    """
}
