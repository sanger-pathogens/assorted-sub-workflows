process MAPPED_READS_TO_FASTQ {
    tag "${meta.ID}"
    label 'cpu_2'
    label 'mem_100M'
    label 'time_30m'
    
    container 'quay.io/biocontainers/samtools:1.22--h96c455f_0'

    input:
    tuple val(meta), path(mapped_reads)

    output:
    tuple val(meta), path("${reads}.gz"), emit: fastq_ch

    script:
    reads = "${meta.ID}.fastq"


    // -f 4 = read unmapped
    // --rf 192 = read is either first in pair or mate in pair
    // pipe to awk brings in line with generate_mags sam_to_fastq.py script
    // (see https://github.com/bxlab/metaWRAP/blob/c4a23f0fff872abe3cafb63767dc592a8361331b/bin/metawrap-scripts/sam_to_fastq.py)
    """
    samtools view \\
        -f 4 \\
        --rf 192 \\
        -@ ${task.cpus} \\
        ${mapped_reads} | \\
    awk -F \$'\t' 'NF >= 11 {print ">"\$1; print \$10; print "+"; print \$11}' \\
    > ${reads}

    gzip ${reads}
    """
}
