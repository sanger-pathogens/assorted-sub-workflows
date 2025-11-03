process SEQKIT {
    tag "${meta.ID}"
    label "cpu_1"
    label "mem_1"
    label "time_30m"

    container  'quay.io/biocontainers/seqkit:2.10.0--h9ee0642_0'

    publishDir mode: 'copy', path: "${params.outdir}/seqkit/"

    input:
    tuple val(meta), path(fastas, stageAs: "fastas/*")

    output:
    tuple val(meta), path(seqkit_out), emit: results

    script:
    seqkit_out = "fasta_out/*"
    """
    mkdir -p fasta_out

    for fasta in fastas/*.${params.fasta_ext}; do
        base=\$(basename "\$fasta")
        seqkit seq "\$fasta" -m ${params.min_contig} > fasta_out/\$base
    done
    """
}
