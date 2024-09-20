process DOWNLOAD_METADATA {
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    container 'quay.io/sangerpathogens/enadownloader:v2.3.2-fb2c2cca'

    input:
    tuple val(meta), path(matches)

    output:
    tuple val(meta), path("metadata.tsv"), emit: metadata_tsv

    script:
    """
    enadownloader --input ${matches} --type sample --write-metadata
    """
}
