// MEGAHIT on the raw R1/R2 pair, for MEGAHIT-only assembly (--metaspades false): the assembly fallback when
// metaSPAdes (and --fastspades) run out of memory/time. A separate process on purpose: MEGAHIT (megahit.nf) takes the
// single file of reads metaSPAdes left unmapped (-r) in the default hybrid mode, and editing its script would change
// its Nextflow cache hash and force every -resume to re-run MEGAHIT and everything after it.
process MEGAHIT_PE {
    tag "${meta.ID}"
    label params.megahit_deterministic ? 'cpu_1' : 'cpu_8'
    label 'mem_8'
    label 'time_12'

    container 'quay.io/biocontainers/megahit:1.2.9--h5ca1c30_6'

    input:
    tuple val(meta), path(first_read), path(second_read)   // same (meta, R1, R2) shape the QC step emits

    output:
    tuple val(meta), path("${meta.ID}_contigs.fasta"), emit: contigs

    script:
    def contigs = "megahit/final.contigs.fa"
    """
    megahit -1 ${first_read} -2 ${second_read} \\
	        -o megahit \\
	        -t ${task.cpus} \\
            -m ${task.memory.toBytes()}

    mv ${contigs} ${meta.ID}_contigs.fasta
    """
}
