process CONVERT_TO_FASTQ {
    label 'cpu_2'
    label 'mem_250M'
    label 'time_30m'

    publishDir path: "${params.outdir}/fastqs/", enabled: params.save_fastqs, mode: 'copy', overwrite: true
    
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1'

    input:
    tuple val(meta), path(reads_bam)

    output:
    tuple val(meta), path(fastq_output),  emit: reads_fastq

    script:
    //fall back to barcode_kit_barcode if no ID
    fastq_output = "${meta.ID ?: "${meta.barcode_kit}_B${meta.barcode}"}.fastq.gz"
    """
    samtools fastq -@ ${task.cpus} -0 ${fastq_output} ${reads_bam}
    """
}

process PUBLISH_BAMS {
    label 'cpu_1'
    label 'mem_1'
    label 'time_30m'

    executor = 'local'

    publishDir path: "${params.outdir}/bams/", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(reads_bam)

    output:
    tuple val(meta), path(bam_output),  emit: reads_bam

    script:
    bam_output = "${meta.ID ?: "${meta.barcode_kit}_B${meta.barcode}"}.bam"
    """
    mv ${reads_bam} ${bam_output}
    """
}
