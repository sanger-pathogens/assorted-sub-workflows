process BUNDLE_FASTAS {
    tag "${meta.ID}"
    input:
    tuple val(meta), path(fasta_list)

    output:
    tuple val(meta), path("fasta_dir")

    script:
    """
    mkdir -p fasta_dir
    cp ${fasta_list.join(' ')} fasta_dir/
    """
}
