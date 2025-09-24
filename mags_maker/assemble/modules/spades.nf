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
def scaffolds = "metaspades/scaffolds.fasta"

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
        -o metaspades \\
        -1 ${first_read} \\
        -2 ${second_read} \\
        \${phred_flag}

mv ${scaffolds} ${meta.ID}_scaffolds.fasta
"""
}