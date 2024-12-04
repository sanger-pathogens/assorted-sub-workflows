process BAKTA {
    tag "${meta.ID}"
    label "cpu_4"
    label "mem_16"
    label "time_1"

    publishDir mode: 'copy', pattern: "${amended_id}.gff3", path: "${params.outdir}/gffs/"

    container 'quay.io/biocontainers/bakta:1.9.4--pyhdfd78af_0'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${amended_id}.embl")             , emit: embl
    tuple val(meta), path("${amended_id}.faa")              , emit: faa
    tuple val(meta), path("${amended_id}.ffn")              , emit: ffn
    tuple val(meta), path("${amended_id}.fna")              , emit: fna
    tuple val(meta), path("${amended_id}.gbff")             , emit: gbk
    tuple val(meta), path("${amended_id}.gff3")             , emit: gff
    tuple val(meta), path("${amended_id}.hypotheticals.tsv"), emit: hypotheticals_tsv
    tuple val(meta), path("${amended_id}.hypotheticals.faa"), emit: hypotheticals_faa
    tuple val(meta), path("${amended_id}.tsv")              , emit: tsv
    tuple val(meta), path("${amended_id}.txt")              , emit: txt
    tuple val(meta), path("${amended_id}.svg")              , emit: svg
    tuple val(meta), path("${amended_id}.png")              , emit: png

    script:
    amended_id = "${meta.ID}".replaceAll(/[^\w.-]/, '_')
    gff = "${meta.ID}.gff3"
    """
    bakta \\
        ${fasta} \\
        ${params.bakta_args} \\
        --threads ${task.cpus} \\
        --prefix ${amended_id} \\
        --locus-tag ${amended_id} \\
        --db ${params.bakta_db} \\
        --keep-contig-headers

    # Remove non-ASCII characters from GFF
    sed -i "s/[‘’]/'/g" "${gff}"
    """
}
