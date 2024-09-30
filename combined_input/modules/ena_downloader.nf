process DOWNLOAD_METADATA {
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    container 'quay.io/sangerpathogens/enadownloader:v2.3.2-fb2c2cca-bookworm'

    if (params.publish_metadata) {
        publishDir "${params.outdir}/metadata", mode: 'copy', overwrite: true
    }

    input:
    tuple val(meta), path(accessions)

    output:
    tuple val(meta), path("metadata.tsv"), emit: metadata_tsv

    script:
    """
    enadownloader --input ${accessions} --type ${params.accession_type} --write-metadata
    """
}
