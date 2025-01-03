process DOWNLOAD_METADATA {
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'
    label 'no_retry'

    container 'quay.io/sangerpathogens/enadownloader:v2.3.3-903be379'

    publishDir "${params.outdir}/metadata", mode: 'copy', overwrite: true, enabled: params.publish_metadata

    input:
    tuple val(meta), path(accessions)

    output:
    tuple val(meta), path("metadata.tsv"), emit: metadata_tsv

    script:
    """
    enadownloader --input ${accessions} --type ${params.accession_type} --write-metadata
    """
}

process DOWNLOAD_FASTQS {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_30m'
    maxForks 10

    publishDir "${params.outdir}/${meta.ID}/fastqs", mode: 'copy', overwrite: true

    input:
    tuple val(meta), val(fastq_path_1), val(fastq_path_2)

    output:
    tuple val(meta), path(read_1), path(read_2), emit: fastqs

    script:
    read_1 = "${meta.ID}_1.fastq.gz"
    read_2 = "${meta.ID}_2.fastq.gz"
    """
    wget --progress=dot:giga \\
            -O ${read_1} \\
            ${fastq_path_1}

    wget --progress=dot:giga \\
            -O ${read_2} \\
            ${fastq_path_2}
    """
}
