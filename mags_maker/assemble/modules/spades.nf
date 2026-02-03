process METASPADES {
    tag "${meta.ID}"
    label 'cpu_8'
    label 'mem_32'
    label 'time_12'

    container 'quay.io/biocontainers/spades:3.15.5--h95f258a_1'

    input:
    tuple val(meta), path(first_read), path(second_read)

    output:
    tuple val(meta), path("${meta.ID}_scaffolds.fasta"), emit: scaffolds

script:
output_folder = "metaspades"
spades_log="${output_folder}/spades.log"
scaffolds = "${output_folder}/scaffolds.fasta"

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

metaspades.py ${params.fastspades ? "--only-assembler" : ""} \\
        --tmp-dir tmp \\
        -t ${task.cpus} \\
        -m ${task.memory.toGiga()} \\
        -o ${output_folder} \\
        -1 ${first_read} \\
        -2 ${second_read} \\
        \${phred_flag}
spades_status=\${?}

if [ \${spades_status} -eq 0 ]; then
    echo "SPAdes completed successfully (exit code \${spades_status})" >&2
    echo "Catching known warnings from spades.log, as there may be issues with the assembly. Known warnings will appear below..." >&2

    ## empty output contigs.fasta file and no scaffold file, often meaning low read input - exit 7
    grep '======= SPAdes pipeline finished WITH WARNINGS!' ${spades_log} 1>&2 \\
        && grep ' * Assembled contigs are in .\\+contigs.fasta' ${spades_log} 1>&2 \\
        && [ ! -s contigs.fasta ] \\
        && exit 7
    
    # NB: in case scaffolds.fasta file is missing not due to the above, nextflow will error out as expecting it as output file
    mv ${scaffolds} ${meta.ID}_scaffolds.fasta

else
    echo "SPAdes failed with exit code \${spades_status}" >&2
    echo "Checking to see if this is a known issue. Any known errors that are found will appear below..." >&2

    ## empty input read file - exit 7
    grep '== Error ==  file is empty' ${spades_log} \\
        && exit 7

    ## "mimalloc: error: unable to allocate OS memory"
    # appears to be exit code 12 in the logs, but spades exits with code 250
    # fixed when sufficient memory is requested, so exit 130 to enable retry strategy
    grep 'mimalloc: error: unable to allocate OS memory' ${spades_log} >&2 \\
        && echo "SPAdes failed due to insufficient memory. Process will be retried with more memory." >&2 \\
        && exit 130

    # sometimes spades catches a memory issues (or at least does not throw the mimalloc error above)
    grep 'mmap(2) failed. Reason: Cannot allocate memory' ${spades_log} >&2 \\
        && echo "SPAdes failed due to insufficient memory. Process will be retried with more memory." >&2 \\
        && exit 130

    ## segmentation fault, possibly due to farm environment and spades not being compiled against the machine/in the singularity container - exit 3
    # This error may not be informative
    grep '== Error ==  system call for:.\\+/usr/local/bin/spades-hammer.\\+finished abnormally' ${spades_log} 1>&2 \\
        && exit 3
fi

# if not caught known exception, process should not have exited yet - do it now with stored metaspades exit status
exit \${spades_status}
"""
}
