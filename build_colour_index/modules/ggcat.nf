process GGCAT {
    tag "${meta.ID}"
    label 'cpu_32'
    label 'mem_64'
    label 'time_queue_from_small'

    container "quay.io/biocontainers/ggcat:2.2.0--hf1b6044_0"

    publishDir mode: 'copy', path: "${params.outdir}/ggcat/${meta.ID}/"

    input:
    tuple val(meta), path(file_colors_input)

    output:
    tuple val(meta), path(unitigs_fna), emit: unitigs

    script:
    unitigs_fna = "unitigs-k${params.kmer_size}.fna"
    """
    ggcat build \\
        -l ${file_colors_input} \\
        -o ${unitigs_fna} \\
        -s 1 \\
        -k ${params.kmer_size} \\
        -t ggcat_temp \\
        -j ${task.cpus} \\
        -m ${task.memory.toGiga()} \\
        -p
    """
}
