process BOWTIE2 {
    label 'cpu_4'
    label 'mem_8'
    label 'time_1'

    conda 'bioconda::bowtie2=2.5.1'
    container 'quay.io/biocontainers/bowtie2:2.5.1--py310h8d7afc0_0'

    input:
    tuple val(meta), path(reads_1), path(reads_2), path(reference), path(bt2_index_files)

    output:
    tuple val(meta), path("${mapped_reads}"),  emit: mapped_reads

    script:
    mapped_reads = "${meta.ID}.sam"
    """
    # glob pattern to ensure correct bt index name
    bt_index=\$(ls *.bt2* | head -1 | awk -F ".1.bt2" '{ print \$1 }')
    bowtie2 -x \${bt_index} \
            -1 ${reads_1} -2 ${reads_2} \
            -S ${mapped_reads} \
            -p ${task.cpus}
    """
}

process BOWTIE2_INDEX {
    label 'cpu_4'
    label 'mem_8'
    label 'time_30m'

    publishDir "${params.outdir}/bowtie2", mode: 'copy', overwrite: true

    conda 'bioconda::bowtie2=2.5.1'
    container 'quay.io/biocontainers/bowtie2:2.5.1--py310h8d7afc0_0'

    input:
    // ref_key is the original reference path as a string; it is passed as a val so
    // it is not staged and survives unchanged, giving downstream joins a stable key
    // (path(reference) below is re-emitted as a work directory path, so it cannot be one)
    tuple val(ref_key), path(reference)

    output:
    tuple val(ref_key), path(reference), path("${reference.baseName}*.bt2"),  emit: bt2_index

    script:
    ref_basename = "${reference.baseName}"
    """
    bowtie2-build ${reference} ${ref_basename}
    """
}
