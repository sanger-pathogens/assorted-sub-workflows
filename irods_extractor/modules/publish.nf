process PUBLISH_FASTQ {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_12'

    publishDir path: { if ("${params.save_method}" == "nested") "${params.outdir}/${meta.ID}/${params.raw_reads_output_dir}/" else "${params.outdir}/${params.raw_reads_output_dir}/" }, enabled: params.save_fastqs, mode: 'copy', overwrite: true, pattern: "${final_name}"

    input:
    tuple val(meta), path(reads_dir)

    output:
    tuple val(meta), path(final_name), emit: path_channel

    script:
    final_name="${meta.ID}.fastq.gz"
    """
    cat ${reads_dir}/*fastq.gz > ${final_name}
    """
}


process PUBLISH_UNBASECALLED {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_12'

    // Only nested publishing for unbasecalled squiggles
    publishDir path: "${params.outdir}/raw_unbasecalled/${meta.ID}/", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(unbasecalled_dir)

    output:
    tuple val(meta), path("*.{pod5,fast5}"), emit: path_channel

    script:
    """
    ln -s \$(readlink -f ${unbasecalled_dir})/*.{pod5,fast5} .
    """
}