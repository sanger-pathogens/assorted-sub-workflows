a
process CONVERT_TO_FASTQ {
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'

    publishDir path: "${params.outdir}/fastqs/", enabled: params.save_fastqs, mode: 'copy', overwrite: true, pattern: "*.fastq.gz"
    
    conda 'bioconda::samtools=1.19'
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1'

    input:
    tuple val(meta), path(reads_bam)

    output:
    tuple val(meta), path(fastq_output),  emit: reads_fastq

    script:
    //will likely need thinking here as if we do other methods of sequencing and use this pipeline need the correct out flags
    fastq_output = "${meta.ID}.fastq.gz"
    """
    samtools fastq -@ ${task.cpus} -0 ${fastq_output} ${reads_bam}
    """
}
a
process MERGE_BAMS_FOR_SUMMARY {
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'
    
    conda 'bioconda::samtools=1.19'
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1'

    input:
    path("*.bam")

    output:
    path(combined_bam),  emit: summary_bam

    script:
    combined_bam = "merged.bam"
    """
    samtools merge -@ ${task.cpus} -o ${combined_bam} *.bam
    """
}
a
process MANAGE_DUPLICATES_FROM_BAMS {
    label 'cpu_2'
    label 'mem_4'
    label 'time_1'

    conda "bioconda::samtools=1.19"
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1"

    input:
    tuple val(meta), path(bam), path(duplicates_list)
    val(mode)

    output:
    tuple val(meta), path(final_bam), emit: bam

    script:
    final_bam = "${bam.simpleName}_clean.bam"
    if (mode == "keep")
        """
        samtools view -N ${duplicates_list} -o ${final_bam} ${bam}
        """
    else if (mode == "remove")
        """
        samtools view -N ^${duplicates_list} -o ${final_bam} ${bam}
        """
}