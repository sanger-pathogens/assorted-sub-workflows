process SUBSET_GTDB {
    label 'cpu_1'
    label 'mem_1'
    label 'time_queue_from_normal'

    input:
    tuple val(meta), path(sourmash_genomes)

    output:
    tuple val(meta), path("subset_gtdb_ref.fasta"), emit: subset_genome
    script:
    """
    while read genome
    do
      zcat ${params.genome_dir}/\${genome}_genomic.fna.gz >> subset_gtdb_ref.fasta
    done < ${sourmash_genomes}
    """
}