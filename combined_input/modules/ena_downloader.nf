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
