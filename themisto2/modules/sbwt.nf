process SBWT_BUILD {
    tag "${meta.ID}"
    label 'cpu_16'
    label 'mem_32'
    label 'time_queue_from_small'

    // Local .sif, not a registry pull -- the registry path (gitlab-registry.internal.sanger.ac.uk/
    // sanger-pathogens/farm_installs/farm-etc/sbwt-rs-cli/0.4.2-f93d92c) was never actually
    // published (404s on tag `latest` under both /apptainer and /docker). This .sif ships
    // alongside the verified-fixed module install instead -- see [[sbwt_setdiff_bugfix_module]]
    // in project memory: 0.4.2-f93d92c (this one) is the module confirmed to contain the actual
    // set-diff fix (2026-08-04); the older 0.4.2-7c5fcc0 looked fixed but was reconfirmed to
    // still corrupt output on 2026-07-31 -- don't swap back to that one.
    container "/data/pam/installs/packages/sbwt-rs-cli/bug_fix_setdiff_commit_f93d92_2026.08.04.13.38.59/sbwt-rs-cli-0.4.2-f93d92c/image/sbwt-rs-cli_bug_fix_setdiff_commit_f93d92_2026.08.04.13.38.59.sif"

    publishDir mode: 'copy', path: "${params.outdir}/sbwt/${meta.ID}/"

    input:
    tuple val(meta), path(unitigs_fna)

    output:
    tuple val(meta), path(sbwt_index), path(lcs_index), emit: index

    script:
    sbwt_index   = "unitigs-k${params.kmer_size}.sbwt"
    lcs_index    = "unitigs-k${params.kmer_size}.lcs"
    // Per-meta.ID subpath -- see ggcat.nf's temp_dir for why (shared params.temp_dir
    // path would otherwise collide across concurrent species/lineage/candidate tasks).
    def temp_dir = params.temp_dir ? "${params.temp_dir}/sbwt/${meta.ID}" : "sbwt_temp"
    // Cap at 95% of the allocated memory -- LSF can kill the job for exceeding its
    // limit even if sbwt asks for exactly what was granted (scheduler overhead eats
    // a little of that headroom), so leave a safety margin.
    // TODO: if the fixed 'mem_32' label ever proves too small for real data (actual
    // OOM-kill, not speculation), see kraken2bracken.config's estimate_kraken_mem()
    // for a pattern that sizes the request from real input + escalates on retry.
    def mem_gb = Math.floor(task.memory.toGiga() * 0.95) as int

    // -l/--build-lcs is required -- without it `sbwt build` only writes the .sbwt file,
    // not the .lcs array, but this process's `output:` above declares lcs_index as
    // mandatory (caused a "Missing output file(s)" failure before this was added).
    """
    mkdir -p ${temp_dir}
    sbwt build \\
        -i ${unitigs_fna} \\
        -o unitigs-k${params.kmer_size} \\
        -r \\
        -l \\
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

    // Same .sif as SBWT_BUILD above -- see that process's comment for why.
    container "/data/pam/installs/packages/sbwt-rs-cli/bug_fix_setdiff_commit_f93d92_2026.08.04.13.38.59/sbwt-rs-cli-0.4.2-f93d92c/image/sbwt-rs-cli_bug_fix_setdiff_commit_f93d92_2026.08.04.13.38.59.sif"

    input:
// Generic -- never touches .lcs (sbwt check only ever takes -i), so this is
    // reused via include-aliasing after both SBWT_BUILD and SBWT_DIFFERENCE.
    // build_color_index.nf re-pairs the checked .sbwt with its .lcs afterward for
    // the Themisto2 build step, since that pairing is plumbing this process
    // doesn't need. 
    tuple val(meta), path(sbwt_index)

    output:
    tuple val(meta), path(sbwt_index), emit: index

    script:
    """
    sbwt check -i ${sbwt_index}
    """
}

process SBWT_DIFFERENCE {
    tag "${meta.ID}"
    label 'cpu_16'
    label 'mem_32' // provisional, sized like SBWT_BUILD -- bg_excl vs ATB needs far more (900GB, hugemem), give it its own withName: override rather than bumping this
    label 'time_queue_from_normal'

    // Same .sif as SBWT_BUILD/SBWT_CHECK -- always pair with SBWT_CHECK on the output (setdiff_filter.nf already does)
    container "/data/pam/installs/packages/sbwt-rs-cli/bug_fix_setdiff_commit_f93d92_2026.08.04.13.38.59/sbwt-rs-cli-0.4.2-f93d92c/image/sbwt-rs-cli_bug_fix_setdiff_commit_f93d92_2026.08.04.13.38.59.sif"

    input:
    tuple val(meta), path(sbwt_a), path(sbwt_b) // sbwt_a - sbwt_b

    output:
    tuple val(meta), path(diff_index), emit: index

    script:
    diff_index = "${meta.ID}.sbwt"
    def low_ram_flag = params.sbwt_diff_low_ram ? "--low-ram" : "" // peak-RAM optimization only, not a correctness fix
    """
    sbwt difference ${sbwt_a} ${sbwt_b} -o ${diff_index} -t ${task.cpus} -v ${low_ram_flag}
    """
}