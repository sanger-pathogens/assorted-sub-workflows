process GGCAT {
    tag "${meta.ID}"
    label 'cpu_32'
    label "mem_8"
    label 'time_queue_from_normal'

    // Only request /tmp space if /tmp is being used (assumes TMPDIR is not set)
    if (!params.temp_dir || params.temp_dir.startsWith("/tmp")) {
        label 'request_temp'
    }

    // scratch used for fast node-local temp storage
    scratch true

    container "quay.io/biocontainers/ggcat:2.2.0--hf1b6044_0"

    // meta.stage (only set for the candidate rebuild) keeps this from colliding with
    // that same lineage's own GGCAT_GROUP output at the same meta.ID -- lineage_index
    // and candidate_index deliberately share meta.ID (the lineage name), required for
    // SET_DIFF_CALCULATIONS' join, so meta.ID alone can't disambiguate them on disk.
    publishDir mode: 'copy', path: "${params.outdir}/ggcat/${meta.stage ? "${meta.stage}/" : ''}${meta.ID}/"

    input:
    tuple val(meta), path(file_colors_input)

    output:
    tuple val(meta), path(unitigs_fna), emit: unitigs

    script:
    // Per-meta.ID subpath -- params.temp_dir is a shared path across every task
    // (species, every lineage, every candidate rebuild), so without this, concurrent
    // GGCAT tasks (e.g. multiple --target_groups lineages building at once) would
    // collide on the same working directory. The no-temp_dir fallback below doesn't
    // need this: it's already relative to each task's own isolated Nextflow work dir.
    def temp_dir = params.temp_dir ? "${params.temp_dir}/ggcat/${meta.ID}" : "ggcat_temp"
    unitigs_fna = "unitigs-k${params.color_index_kmer_size}.fna"
    // Cap at 95% of the allocated memory -- LSF can kill the job for exceeding its
    // limit even if ggcat asks for exactly what was granted (scheduler overhead eats
    // a little of that headroom), so leave a safety margin.
    // TODO: mem_8 escalates to 16/32.GB on retry -- if that's ever not enough,
    // see kraken2bracken.config's estimate_kraken_mem() for input-based sizing.
    def mem_gb = Math.floor(task.memory.toGiga() * 0.95) as int
    """
    mkdir -p ${temp_dir}
    ggcat build \\
        -l ${file_colors_input} \\
        -o ${unitigs_fna} \\
        -s 1 \\
        -k ${params.color_index_kmer_size} \\
        -t ${temp_dir} \\
        -j ${task.cpus} \\
        -m ${mem_gb} \\
        -p
    """
}