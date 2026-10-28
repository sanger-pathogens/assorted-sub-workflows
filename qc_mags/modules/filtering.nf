process FILTER_FASTAS {
    tag "${meta.ID}"
    label "cpu_1"
    label "mem_1"
    label "time_30m"

    container  'quay.io/biocontainers/seqkit:2.10.0--h9ee0642_0'

    publishDir mode: 'copy', path: "${params.outdir}/pass/fastas"

    input:
    tuple val(meta), val(filenames), path(fastas, stageAs: "fastas/*")

    output:
    tuple val(meta), path("*.${fasta_ext}"), emit: filtered_fastas

    script:
    fasta_ext = params.fasta_ext.replaceAll(/^\./, '')
    """
    for file in ${filenames.join(" ")}; do
        ln -s fastas/\${file} .
    done
    """
}

process PUBLISH_RESULTS {
    label "cpu_1"
    label "mem_100M"
    label "time_1"

    publishDir mode: 'copy', path: "${output_dir}"

    input:
    tuple val(meta), path(inputs)
    val(output_dir)

    output:
    tuple val(meta), path(inputs)

    script:
    """
    echo ${output_dir} > published.txt
    """
}
