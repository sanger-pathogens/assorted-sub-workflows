process MEGAHIT {
    tag "${meta.ID}"
    label params.megahit_deterministic ? 'cpu_1' : 'cpu_8'
    label 'mem_8'
    label 'time_12'

    container 'quay.io/biocontainers/megahit:1.2.9--h5ca1c30_6'

    input:
    tuple val(meta), path(unmapped_reads)

    output:
    tuple val(meta), path("${meta.ID}_contigs.fasta"), emit: contigs

    script:
    def contigs = "megahit/final.contigs.fa"
    """
    megahit -r ${unmapped_reads} \\
	        -o megahit \\
	        -t ${task.cpus} \\
            -m ${task.memory.toBytes()}

    mv ${contigs} ${meta.ID}_contigs.fasta
    """
}