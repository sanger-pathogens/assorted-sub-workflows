process SAM_TO_FASTQ {
    tag "${meta.ID}"
    label 'cpu_2'
    label 'mem_100M'
    label 'time_30m'
    
    container 'quay.io/biocontainers/samtools:1.22--h96c455f_0'

    input:
    tuple val(meta), path(sam)

    output:
    tuple val(meta), path("${reads}.gz"), emit: fastq_ch

    script:
    reads = "${meta.ID}.fastq"


    //-f 4 = read unmapped
    //--rf 192 = read is either first in pair or mate in pair
    """
    samtools view -b \\
        -f 4 \\
        --rf 192 \\
        -@ ${task.cpus} \\
        ${sam} | \\
    samtools fastq -N \\
        -o ${reads} \\
        -@ ${task.cpus} \\
        -

    gzip ${reads}
    """
}