process BOWTIE2 {
    label 'cpu_4'
    label 'mem_8'
    label 'time_12'

    publishDir "${params.outdir}/${meta.id}/bowtie2", mode: 'copy', overwrite: true

    container 'quay.io/biocontainers/bowtie2:2.5.1--py310h8d7afc0_0'

    input:
    tuple val(meta), file(reads)
    path(bt2_index_files)

    output:
    tuple val(meta), path("${mapped_reads}"),  emit: mapped_reads

    script:
    mapped_reads = "${meta.id}.sam"
    """
    # glob pattern to ensure correct bt index name
    bt_index=\$(ls *.bt2* | head -1 | awk -F ".1.bt2" '{ print \$1 }')
    bowtie2 -x \${bt_index} \
            -1 ${reads[0]} -2 ${reads[1]} \
            -S ${mapped_reads} \
            -p ${task.cpus}
    """
}

process BOWTIE2_INDEX {
    label 'cpu_4'
    label 'mem_8'
    label 'time_12'

    publishDir "${params.outdir}/bowtie2", mode: 'copy', overwrite: true

    container 'quay.io/biocontainers/bowtie2:2.5.1--py310h8d7afc0_0'

    input:
    path(reference)

    output:
    path("${reference.baseName}*.bt2"),  emit: bt2_index

    script:
    ref_basename = "${reference.baseName}"
    """
    bowtie2-build ${reference} ${ref_basename}
    """
}
