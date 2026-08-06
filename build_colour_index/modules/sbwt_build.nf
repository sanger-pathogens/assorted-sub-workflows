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
    def temp_dir = params.sbwt_temp_dir ?: "sbwt_temp"
    
    """
    mkdir -p ${temp_dir}
    sbwt build \\
        -i ${unitigs_fna} \\
        -o unitigs-k${params.kmer_size} \\
        -r \\
        -k ${params.kmer_size} \\
        -m ${task.memory.toGiga()} \\
        -t ${task.cpus} \\
        --temp-dir ${temp_dir}

    sbwt check -i ${sbwt_index}
    """
}
