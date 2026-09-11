process GGCAT {
    tag "${meta.ID}"
    label 'cpu_32'
    label 'time_queue_from_normal'

    // Set directly here, not via a config withName override: GGCAT_SPECIES scales with
    // genome count (measured ~17 GB at ~42k S. pneumoniae genomes; 24 GB baseline leaves
    // headroom), GGCAT_CANDIDATE's input is a handful of candidate marker unitigs --
    // tiny regardless of species size. meta.stage distinguishes the two aliases (unset for
    // species-wide, 'candidate' for the rebuild -- see build_color_index.nf/marker_filtering.nf).
    memory = { meta.stage == 'candidate' ? (8.GB * task.attempt) : (24.GB * task.attempt) }

    // request /tmp only if /tmp is actually the temp dir (assumes TMPDIR unset)
    if (!params.temp_dir || params.temp_dir.startsWith("/tmp")) {
        label 'request_temp'
    }
    scratch true

    container "quay.io/biocontainers/ggcat:2.2.0--hf1b6044_0"

    // meta.stage disambiguates the candidate rebuild (GGCAT_CANDIDATE) from the
    // species-wide GGCAT_SPECIES output ('ggcat/candidate_<lineage>/' vs 'ggcat/<species>/').
    publishDir mode: 'copy', path: "${params.outdir}/ggcat/${meta.stage ? "${meta.stage}_" : ''}${meta.ID}/", enabled: params.publish_intermediate

    input:
    tuple val(meta), path(file_colors_input)

    output:
    tuple val(meta), path(unitigs_fna), emit: unitigs

    script:
    // Per-meta.ID subpath so concurrent tasks don't collide on a shared temp_dir.
    def temp_dir = params.temp_dir ? "${params.temp_dir}/ggcat/${meta.ID}" : "ggcat_temp"
    unitigs_fna = "unitigs-k${params.color_index_kmer_size}.fna"
    // Cap -m at 95% of the allocation -- LSF kills a job that hits its exact limit.
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
        -p \\
        -e
    """
}