process CONVERT_TO_BAM {
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'
    
    conda 'bioconda::samtools=1.17'
    container 'quay.io/biocontainers/samtools:1.17--hd87286a_2'

    input:
    tuple val(meta), path(mapped_reads)

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
    container 'quay.io/biocontainers/samtools:1.17--hd87286a_2'

    input:
    tuple val(meta), path(mapped_reads_bam)

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
    label 'cpu_1'
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

process INDEX_BAM {
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'

    // this publish statement might duplicate output from deeptools bigwig, 
    // but given the small size of index file and the fact that it's most handy when located in same folder as .bam or .bw, it is best publishing it twice
    publishDir "${params.outdir}/${meta.ID}/samtools_sort", enabled: params.keep_sorted_bam, mode: 'copy', overwrite: true

    conda 'bioconda::samtools=1.17'
    container 'quay.io/biocontainers/samtools:1.17--hd87286a_2'

    input:
    tuple val(meta), path(mapped_reads_bam)

    output:
    tuple val(meta), path(mapped_reads_bam), path(mapped_reads_bai),  emit: sorted_indexed_bam

    script:
    mapped_reads_prefix = mapped_reads_bam.simpleName
    mapped_reads_bai = "${mapped_reads_prefix}.bai"
    """
    samtools index -@ ${task.cpus} "${mapped_reads_bam}" "${mapped_reads_bai}"
    """
}

process SAMTOOLS_STATS {
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'

    publishDir "${params.outdir}/${meta.ID}/samtools_stats", enabled: params.samtools_stats, mode: 'copy', overwrite: true

    conda 'bioconda::samtools=1.17'
    container 'quay.io/biocontainers/samtools:1.17--hd87286a_2'

    input:
    tuple val(meta), path(mapped_reads_bam), path(mapped_reads_bai)

    output:
    tuple path(stats_file), path(flagstats_file),  emit: stats_ch

    script:
    stats_file = "${meta.ID}.stats"
    flagstats_file = "${meta.ID}.flagstats"
    """
    samtools stats -@ ${task.cpus} "${mapped_reads_bam}" > "${stats_file}"
    samtools flagstats -@ ${task.cpus} "${mapped_reads_bam}" > "${flagstats_file}"
    """
}

