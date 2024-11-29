process BAKTA {
    tag "${meta.ID}"
    label "cpu_4"
    label "mem_16"
    label "time_1"

    publishDir mode: 'copy', pattern: "${meta.ID}.gff3", path: "${params.outdir}/gffs/"

    container 'quay.io/biocontainers/bakta:1.9.4--pyhdfd78af_0'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.ID}.embl")             , emit: embl
    tuple val(meta), path("${meta.ID}.faa")              , emit: faa
    tuple val(meta), path("${meta.ID}.ffn")              , emit: ffn
    tuple val(meta), path("${meta.ID}.fna")              , emit: fna
    tuple val(meta), path("${meta.ID}.gbff")             , emit: gbk
    tuple val(meta), path("${meta.ID}.gff3")             , emit: gff
    tuple val(meta), path("${meta.ID}.hypotheticals.tsv"), emit: hypotheticals_tsv
    tuple val(meta), path("${meta.ID}.hypotheticals.faa"), emit: hypotheticals_faa
    tuple val(meta), path("${meta.ID}.tsv")              , emit: tsv
    tuple val(meta), path("${meta.ID}.txt")              , emit: txt
    tuple val(meta), path("${meta.ID}.svg")              , emit: svg
    tuple val(meta), path("${meta.ID}.png")              , emit: png

    script:
    """
    bakta \\
        ${fasta} \\
        ${params.bakta_args} \\
        --threads ${task.cpus} \\
        --prefix ${meta.ID} \\
        --locus-tag ${meta.bakta_id} \\
        --db ${params.bakta_db} \\
        --keep-contig-headers
    """
}