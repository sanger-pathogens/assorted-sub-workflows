process SEQKIT {
    tag "${meta.ID}"
    label "cpu_1"
    label "mem_1"
    label "time_30m"

    container  'quay.io/biocontainers/seqkit:2.10.0--h9ee0642_0'

    publishDir mode: 'copy', path: "${params.outdir}/seqkit/"

    input:
    tuple val(group_key), val(meta), path(fasta)

    output:
    tuple val(group_key), val(meta), path(finalName), emit: results

    script:
    finalName = "cleaned_${fasta.getBaseName()}" 
    """
    seqkit seq ${fasta} -m ${params.min_contig_length} > ${finalName}
    """
}