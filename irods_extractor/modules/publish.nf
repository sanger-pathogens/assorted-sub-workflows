process PUBLISH_FASTQ {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_30m'

    publishDir path: { if ("${params.save_method}" == "nested") "${params.outdir}/${meta.ID}/${params.raw_reads_prefix}fastqs/" else "${params.outdir}/${params.raw_reads_prefix}fastqs/" }, enabled: params.save_fastqs, mode: 'copy', overwrite: true, pattern: "${final_name}"

    input:
    val(meta)

    output:
    tuple val(meta), path(final_name), emit: path_channel

    script:
    //fix for pod5?
    final_name="${meta.ID}.fastq.gz"
    """
    cat ${meta.local_path}/*fastq.gz > ${final_name}
    """
}