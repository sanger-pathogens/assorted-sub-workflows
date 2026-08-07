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
    // Cap at 95% of the allocated memory -- LSF can kill the job for exceeding its
    // limit even if ggcat asks for exactly what was granted (scheduler overhead eats
    // a little of that headroom), so leave a safety margin.
    // TODO: if the fixed 'mem_64' label ever proves too small for real data (actual
    // OOM-kill, not speculation), see kraken2bracken.config's estimate_kraken_mem()
    // for a pattern that sizes the request from real input + escalates on retry.
    def mem_gb = Math.floor(task.memory.toGiga() * 0.95) as int
    """
    ggcat build \\
        -l ${file_colors_input} \\
        -o ${unitigs_fna} \\
        -s 1 \\
        -k ${params.kmer_size} \\
        -t ${params.ggcat_temp_dir} \\
        -j ${task.cpus} \\
        -m ${mem_gb} \\
        -p
    """
}