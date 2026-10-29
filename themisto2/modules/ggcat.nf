process GGCAT {
    tag "${meta.ID}"
    label 'cpu_32'
    label "mem_8"
    label 'time_queue_from_normal'

    // request /tmp only if /tmp is actually the temp dir (assumes TMPDIR unset)
    if (!params.temp_dir || params.temp_dir.startsWith("/tmp")) {
        label 'request_temp'
    }
    scratch true

    container "quay.io/biocontainers/ggcat:2.2.0--hf1b6044_0"

    // meta.stage disambiguates the candidate rebuild (GGCAT_CANDIDATE) from the
    // species-wide GGCAT_SPECIES output when they share a publishDir root.
    publishDir mode: 'copy', path: "${params.outdir}/ggcat/${meta.stage ? "${meta.stage}/" : ''}${meta.ID}/", enabled: params.publish_intermediate

    input:
    tuple val(meta), path(file_colors_input)

    output:
    tuple val(meta), path(unitigs_fna), emit: unitigs

    script:
    // Per-meta.ID subpath so concurrent tasks don't collide on a shared temp_dir.
    def temp_dir = params.temp_dir ? "${params.temp_dir}/ggcat/${meta.ID}" : "ggcat_temp"
    unitigs_fna = "unitigs-k${params.color_index_kmer_size}.fna"
    // Provide memory target 10% lower than LSF job suggested by Jarno
    def mem_gb = (task.memory.toGiga() * 0.9).toInteger()
    """
    mkdir -p ${temp_dir}
    ggcat build \\
        --input-lists ${file_colors_input} \\
        --output-file ${unitigs_fna} \\
        --min-multiplicity 1 \\
        --kmer-length ${params.color_index_kmer_size} \\
        --temp-dir ${temp_dir} \\
        --threads-count ${task.cpus} \\
        --memory ${mem_gb} \\
        --prefer-memory
    """
}