process SPADES_REASSEMBLE {
    tag "${meta.ID}_${fullBinInfo.bin}_${fullBinInfo.level}"
    label 'cpu_8'
    label 'mem_32'
    label 'time_12'

    container 'quay.io/biocontainers/spades:3.15.5--h95f258a_1'

    input:
    tuple val(meta), val(fullBinInfo), path(bin), path(first_read), path(second_read)

    output:
    tuple val(meta), path(final_name), emit: contigs

    script:
    def contigs = "reassembled/contigs.fasta"
    final_name = "${meta.ID}_bin_${fullBinInfo.bin}_${fullBinInfo.level}.fasta"
    """
    spades.py \\
            --tmp-dir tmp \\
            -t ${task.cpus} \\
            -m ${task.memory.toGiga()} \\
            --careful \\
            --untrusted-contigs ${bin} \\
            -o reassembled \\
            -1 ${first_read} \\
            -2 ${second_read} \\
            ${params.lock_phred ? "--phred-offset 33" : ""}

    mv ${contigs} ${final_name}
    """
}