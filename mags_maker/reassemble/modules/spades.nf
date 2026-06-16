process SPADES_REASSEMBLE {
    tag "${meta.ID}_${fullBinInfo.bin}_${fullBinInfo.level}"
    label 'cpu_8'
    label 'mem_32'
    label 'time_12'

    container 'quay.io/biocontainers/spades:3.15.5--h95f258a_1'

    input:
    tuple val(meta), val(fullBinInfo), path(bin), path(first_read), path(second_read)

    output:
    tuple val(meta), path(long_scaffolds), emit: long_scaffolds

    script:
    def scaffolds = "reassembled/scaffolds.fasta"
    final_name = "${meta.ID}_bin_${fullBinInfo.bin}_${fullBinInfo.level}.fasta"
    def long_scaffolds = "long_${final_name.name}"
    def careful = params.careful ? "--careful" : ""

    command = "${projectDir}/assorted-sub-workflows/mags_maker/assemble/bin/rm_short_contigs.py"
    min_contig_length = [params.maxbin2_min_contig, params.concoct_min_contig, params.metabat_min_contig].min()
    """
    # This is done because if the sra-lite format there is no quality information so --phred-offset needs to be set
    # Determine phred flag

    if [[ "${params.lock_phred}" == "true" ]]; then
        phred_flag="--phred-offset 33"
    elif grep -q '?' <(zcat "${first_read}" | head -n 75); then
        phred_flag="--phred-offset 33"
    else
        phred_flag=""
    fi

    spades.py \\
            --tmp-dir tmp \\
            -t ${task.cpus} \\
            -m ${task.memory.toGiga()} \\
            ${careful} \\
            --untrusted-contigs ${bin} \\
            -o reassembled \\
            -1 ${first_read} \\
            -2 ${second_read} \\
            \${phred_flag}

    mv ${scaffolds} ${final_name}

    status=\${?}
    if [ \${status} -gt 0 ] ; then
        # remap exit 12 memory error to 130 to enable retry strategy
        grep 'Cannot allocate memory. Error code: 12' spades.log 
        exit 130
    else
        exit \$status
    fi

    # Remove small contigs
    ${command} ${min_contig_length} ${contigs} > ${long_scaffolds}

    """
}