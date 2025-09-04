process SAM_TO_FASTQ {
    label 'cpu_2'
    label 'mem_100M'
    label 'time_30m'
    
    container 'quay.io/biocontainers/samtools:1.22--h96c455f_0'

    input:
    tuple val(meta), path(sam)

    output:
    tuple val(meta), path("${read1}.gz"), path("${read2}.gz"), emit: fastq_ch

    script:
    read1 = "${meta.ID}_1.fastq"
    read2 = "${meta.ID}_2.fastq"

    //-f 4 = read unmapped
    //--rf 192 = read is either first in pair or mate in pair
    """
    samtools view -b \\
        -f 4 \\
        --rf 192 \\
        -@ ${task.cpus} \\
        ${sam} | \\
    samtools fastq -N \\
        -1 ${read1} \\
        -2 ${read2} \\
        -@ ${task.cpus} \\
        -

    gzip ${read1} ${read2}
    """
}