process CONVERT_TO_BAM {
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'
    
    conda 'bioconda::samtools=1.17'
    container "${ profile.name == 'standard' ? '/software/pathogen/images/samtools-1.17.simg' : 'quay.io/biocontainers/samtools:1.17--hd87286a_2' }"

    input:
    tuple val(meta), file(mapped_reads)

    output:
    tuple val(meta), path("${mapped_reads_bam}"),  emit: mapped_reads_bam

    script:
    mapped_reads_bam = "${meta.ID}.bam"
    """
    samtools view -@ ${task.cpus} \
                  -bS -F4 \
                  -o ${mapped_reads_bam} \
                  ${mapped_reads}
    """
}

process SAMTOOLS_SORT {
    label 'cpu_4'
    label 'mem_8'
    label 'time_12'

    publishDir "${params.outdir}/${meta.ID}/samtools_sort", enabled: params.keep_sorted_bam, mode: 'copy', overwrite: true

    conda 'bioconda::samtools=1.17'
    container "${ profile.name == 'standard' ? '/software/pathogen/images/samtools-1.17.simg' : 'quay.io/biocontainers/samtools:1.17--hd87286a_2' }"

    input:
    tuple val(meta), file(mapped_reads_bam)

    output:
    tuple val(meta), path("${sorted_reads}"),  emit: sorted_reads

    script:
    sorted_reads = "${meta.ID}_sorted.bam"
    """
    samtools sort -@ ${task.cpus} \
                  -o ${sorted_reads} \
                  ${mapped_reads_bam}
    """
}

process INDEX_REF {
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'

    publishDir "${params.outdir}/sorted_ref", mode: 'copy', overwrite: true

    conda 'bioconda::samtools=1.17'
    container 'quay.io/biocontainers/samtools:1.17--hd87286a_2'

    input:
    path(reference)

    output:
    tuple path(reference), path("${faidx}"),  emit: ref_index

    script:
    faidx = "${reference}.fai"
    """
    samtools faidx "${reference}" > "${faidx}"
    """
}
