// All SBWT_* processes use a local .sif, not a registry pull: the bug-fixed
// sbwt-rs-cli 0.4.2-f93d92c (set-diff corruption fix, 2026-08-04) was never
// published to a registry. Do not swap back to 0.4.2-7c5fcc0 -- still corrupts.
// See [[sbwt_setdiff_bugfix_module]]. Tracked for a public image: PAT-3572.

process SBWT_BUILD {
    tag "${meta.ID}"
    label 'cpu_32'
    label 'mem_8'
    label 'time_queue_from_normal'

    if (!params.temp_dir || params.temp_dir.startsWith("/tmp")) {
        label 'request_temp'
    }
    scratch true

    container "/data/pam/installs/packages/sbwt-rs-cli/bug_fix_setdiff_commit_f93d92_2026.08.04.13.38.59/sbwt-rs-cli-0.4.2-f93d92c/image/sbwt-rs-cli_bug_fix_setdiff_commit_f93d92_2026.08.04.13.38.59.sif"

    publishDir mode: 'copy', path: "${params.outdir}/sbwt/${meta.stage ? "${meta.stage}_" : ''}${meta.ID}/", enabled: params.publish_intermediate

    input:
    tuple val(meta), path(unitigs_fna)

    output:
    tuple val(meta), path(sbwt_index), path(lcs_index), emit: index

    script:
    sbwt_index   = "unitigs-k${params.color_index_kmer_size}.sbwt"
    lcs_index    = "unitigs-k${params.color_index_kmer_size}.lcs"
    def temp_dir = params.temp_dir ? "${params.temp_dir}/sbwt/${meta.ID}" : "sbwt_temp"
    // Cap -m at 95% of the allocation -- LSF kills a job that hits its exact limit.
    def mem_gb = Math.floor(task.memory.toGiga() * 0.95) as int
    // -l is required for the .lcs output declared above (else "Missing output file(s)").
    """
    mkdir -p ${temp_dir}
    sbwt build \\
        -i ${unitigs_fna} \\
        -o unitigs-k${params.color_index_kmer_size} \\
        -r \\
        -l \\
        -k ${params.color_index_kmer_size} \\
        -m ${mem_gb} \\
        -t ${task.cpus} \\
        --temp-dir ${temp_dir}
    """
}

process SBWT_CHECK {
    tag "${meta.ID}"
    label 'cpu_4'
    label 'mem_2'
    label 'time_30m'

    container "/data/pam/installs/packages/sbwt-rs-cli/bug_fix_setdiff_commit_f93d92_2026.08.04.13.38.59/sbwt-rs-cli-0.4.2-f93d92c/image/sbwt-rs-cli_bug_fix_setdiff_commit_f93d92_2026.08.04.13.38.59.sif"

    publishDir mode: 'copy', path: {
        meta.ID.startsWith('markers_')
            ? "${params.outdir}/${meta.ID}/"
            : "${params.outdir}/sbwt/${meta.stage ? "${meta.stage}_" : ''}${meta.ID}/"
    }

    input:
    tuple val(meta), path(sbwt_index)

    output:
    tuple val(meta), path(sbwt_index), emit: index

    script:
    """
    sbwt check -i ${sbwt_index} -t ${task.cpus}
    """
}

process SBWT_DIFFERENCE {
    tag "${meta.ID}"
    label 'cpu_16'
    label 'mem_8'
    label 'time_queue_from_normal'

    container "/data/pam/installs/packages/sbwt-rs-cli/bug_fix_setdiff_commit_f93d92_2026.08.04.13.38.59/sbwt-rs-cli-0.4.2-f93d92c/image/sbwt-rs-cli_bug_fix_setdiff_commit_f93d92_2026.08.04.13.38.59.sif"

    input:
    tuple val(meta), path(sbwt_a, stageAs: 'a.sbwt'), path(sbwt_b, stageAs: 'b.sbwt') // a - b

    output:
    tuple val(meta), path(diff_index), emit: index

    script:
    diff_index = "${meta.ID}.sbwt"
    def low_ram_flag = params.sbwt_diff_low_ram ? "--low-ram" : "" // peak-RAM only, not a correctness fix
    """
    sbwt difference a.sbwt b.sbwt -o ${diff_index} -t ${task.cpus} -v ${low_ram_flag}
    """
}

process SBWT_DUMP_UNITIGS {
    tag "${meta.ID}"
    label 'cpu_4'
    label 'mem_2'
    label 'time_30m'

    container "/data/pam/installs/packages/sbwt-rs-cli/bug_fix_setdiff_commit_f93d92_2026.08.04.13.38.59/sbwt-rs-cli-0.4.2-f93d92c/image/sbwt-rs-cli_bug_fix_setdiff_commit_f93d92_2026.08.04.13.38.59.sif"

    publishDir mode: 'copy', path: "${params.outdir}/${meta.ID}/"

    input:
    tuple val(meta), path(sbwt_index)

    output:
    tuple val(meta), path(unitigs_fasta), emit: unitigs

    script:
    unitigs_fasta = "${meta.ID}_unitigs.fasta"
    """
    sbwt dump-unitigs -i ${sbwt_index} -o ${unitigs_fasta} -t ${task.cpus} -v
    """
}
