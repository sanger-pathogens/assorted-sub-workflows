process RETRIEVE_CRAM {
    label 'cpu_1'
    label 'mem_1'
    label 'time_30m'

    //using the singularity image below the ISG/experiemntal/irods/4.3.0 module
    container '/software/isg/private/experimental/irods/4.3.0/lib/4.3.0_135981_feat-inital.sif'

    input:
    tuple val(meta), val(cram_path)

    output:
    tuple val(meta), path(output_file), emit: path_channel

    script:
    output_file = file(cram_path).name
    """
    iget -K ${cram_path}
    """
}

process RETRIEVE_FASTQ {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_30m'

    publishDir path: { if ("${params.save_method}" == "nested") "${params.outdir}/${meta.ID}/${params.raw_reads_prefix}fastqs/" else "${params.outdir}/fastqs/" }, enabled: params.save_fastqs, mode: 'copy', overwrite: true, pattern: "${final_name}"

    //using the singularity image below the ISG/experiemntal/irods/4.3.0 module
    container '/software/isg/private/experimental/irods/4.3.0/lib/4.3.0_135981_feat-inital.sif'

    input:
    tuple val(meta), val(barcode_path)

    output:
    tuple val(meta), path(final_name), emit: path_channel

    script:
    final_name="${meta.ID}.fastq.gz"
    """
    iget -r  ${barcode_path}
    cat barcode*/*fastq.gz > ${final_name}
    """
}