process BOWTIE2 {
    label 'cpu_4'
    label 'mem_8'
    label 'time_12'

    conda 'bioconda::bowtie2=2.5.1'
    container "${ singularity.enabled ? '/software/pathogen/images/bowtie2-2.5.1--py38he00c5e5_2.simg' : 'quay.io/biocontainers/bowtie2:2.5.1--py310h8d7afc0_0' }"

    input:
    tuple val(meta), path(reads_1), path(reads_2)
    path(bt2_index_files)

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
    label 'time_12'

    publishDir "${params.outdir}/bowtie2", mode: 'copy', overwrite: true

    conda 'bioconda::bowtie2=2.5.1'
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
