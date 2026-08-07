process SBWT_BUILD {
    tag "${meta.ID}"
    label 'cpu_16'
    label 'mem_32'
    label 'time_queue_from_small'

    // Draft registry path -- unconfirmed as the right way to ship this, per team check-in.
    container "gitlab-registry.internal.sanger.ac.uk/sanger-pathogens/farm_installs/farm-etc/sbwt-rs-cli/0.4.2-f93d92c/apptainer"

    publishDir mode: 'copy', path: "${params.outdir}/sbwt/${meta.ID}/"

    input:
    tuple val(meta), path(unitigs_fna)

    output:
    tuple val(meta), path(sbwt_index), path(lcs_index), emit: index

    script:
    sbwt_index   = "unitigs-k${params.kmer_size}.sbwt"
    lcs_index    = "unitigs-k${params.kmer_size}.lcs"
    def temp_dir = params.temp_dir ? "${params.temp_dir}/sbwt" : "sbwt_temp"
    // Cap at 95% of the allocated memory -- LSF can kill the job for exceeding its
    // limit even if sbwt asks for exactly what was granted (scheduler overhead eats
    // a little of that headroom), so leave a safety margin.
    // TODO: if the fixed 'mem_32' label ever proves too small for real data (actual
    // OOM-kill, not speculation), see kraken2bracken.config's estimate_kraken_mem()
    // for a pattern that sizes the request from real input + escalates on retry.
    def mem_gb = Math.floor(task.memory.toGiga() * 0.95) as int

    """
    mkdir -p ${temp_dir}
    sbwt build \\
        -i ${unitigs_fna} \\
        -o unitigs-k${params.kmer_size} \\
        -r \\
        -k ${params.kmer_size} \\
        -m ${mem_gb} \\
        -t ${task.cpus} \\
        --temp-dir ${temp_dir}
    """
}

process SBWT_CHECK {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_16' // provisional guess -- much lighter workload than the build, revisit after seeing real usage in test runs
    label 'time_queue_from_small'

    container "gitlab-registry.internal.sanger.ac.uk/sanger-pathogens/farm_installs/farm-etc/sbwt-rs-cli/0.4.2-f93d92c/apptainer"

    input:
    tuple val(meta), path(sbwt_index), path(lcs_index)

    output:
    tuple val(meta), path(sbwt_index), path(lcs_index), emit: index

    script:
    """
    sbwt check -i ${sbwt_index}
    """
}

// TODO: SBWT_DIFFERENCE process goes here -- will be called by the upcoming
// set_diff_calculations.nf subworkflow. Note from experience: sbwt diff jobs have
// taken a lot of memory when calculating -- try the sbwt `--low-ram` option during
// test runs and see if it helps before locking in a memory label for this process.

//process SBWT_DIFFERENCE {
//}