process BAM_COVERAGE {
    label 'cpu_4'
    label 'mem_16'
    label 'time_30m'

    container "quay.io/biocontainers/deeptools:3.5.2--pyhdfd78af_1"
    
    publishDir "${params.outdir}/${meta.ID}/deeptools_bigwigs/", mode: 'copy', overwrite: true, pattern: "${bigwig}"
    publishDir "${params.outdir}/${meta.ID}/deeptools_bigwigs/", mode: 'copy', overwrite: true, pattern: "${bam_index}"

    input:
    tuple val(meta), path(bam_file), path(bam_index)

    output:
    tuple val(meta), path(bam_index), path(bigwig), emit: bigwig_channel

    script:
    bigwig="${meta.ID}.bw"
    """
    bamCoverage -b ${bam_file} -o ${bigwig}
    """
}