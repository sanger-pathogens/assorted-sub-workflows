process SBWT_BUILD {
    tag "${meta.ID}"
    label 'cpu_32'
    label 'mem_8'
    label 'time_queue_from_normal'

    // Only request /tmp space if /tmp is being used (assumes TMPDIR is not set)
    if (!params.temp_dir || params.temp_dir.startsWith("/tmp")) {
        label 'request_temp'
    }

    // scratch used for fast node-local temp storage
    scratch true

    // Local .sif, not a registry pull -- the registry path (gitlab-registry.internal.sanger.ac.uk/
    // sanger-pathogens/farm_installs/farm-etc/sbwt-rs-cli/0.4.2-f93d92c) was never actually
    // published (404s on tag `latest` under both /apptainer and /docker). This .sif ships
    // alongside the verified-fixed module install instead -- see [[sbwt_setdiff_bugfix_module]]
    // in project memory: 0.4.2-f93d92c (this one) is the module confirmed to contain the actual
    // set-diff fix (2026-08-04); the older 0.4.2-7c5fcc0 looked fixed but was reconfirmed to
    // still corrupt output on 2026-07-31 -- don't swap back to that one.
    container "/data/pam/installs/packages/sbwt-rs-cli/bug_fix_setdiff_commit_f93d92_2026.08.04.13.38.59/sbwt-rs-cli-0.4.2-f93d92c/image/sbwt-rs-cli_bug_fix_setdiff_commit_f93d92_2026.08.04.13.38.59.sif"

    // meta.stage disambiguates the candidate rebuild from that same lineage's own
    // SBWT_BUILD_GROUP output -- see ggcat.nf's publishDir for why meta.ID alone can't.
    publishDir mode: 'copy', path: "${params.outdir}/sbwt/${meta.stage ? "${meta.stage}/" : ''}${meta.ID}/"

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
    // TODO: mem_8 escalates to 16/32.GB on retry -- if that's ever not enough,
    // see kraken2bracken.config's estimate_kraken_mem() for input-based sizing.
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
    label 'cpu_4'
    label 'mem_2' // sized from real test-run peak RSS (<14 MB) -- much lighter than the build, as expected
    label 'time_30m'

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
    sbwt check -i ${sbwt_index} -t ${task.cpus}
    """
}

process SBWT_DIFFERENCE {
    tag "${meta.ID}"
    label 'cpu_16'
    // Sized for XLIN_BG/LIN_CAND (<1.3GB peak, real runs). MARKERS/BG_EXCL diff
    // against ATB-scale data instead -- see setdiff_filter.config's withName override.
    label 'mem_8'
    label 'time_queue_from_normal'

    // Same .sif as SBWT_BUILD/SBWT_CHECK -- always pair with SBWT_CHECK on the output (setdiff_filter.nf already does)
    container "/data/pam/installs/packages/sbwt-rs-cli/bug_fix_setdiff_commit_f93d92_2026.08.04.13.38.59/sbwt-rs-cli-0.4.2-f93d92c/image/sbwt-rs-cli_bug_fix_setdiff_commit_f93d92_2026.08.04.13.38.59.sif"

    input:
    // stageAs required: every SBWT_BUILD output is named unitigs-k${kmer_size}.sbwt
    // regardless of species/lineage/group, so sbwt_a and sbwt_b routinely arrive
    // with the same basename -- without distinct stageAs names Nextflow can't
    // stage both into the task dir ("input file name collision").
    tuple val(meta), path(sbwt_a, stageAs: 'a.sbwt'), path(sbwt_b, stageAs: 'b.sbwt') // sbwt_a - sbwt_b

    output:
    tuple val(meta), path(diff_index), emit: index

    script:
    diff_index = "${meta.ID}.sbwt"
    def low_ram_flag = params.sbwt_diff_low_ram ? "--low-ram" : "" // peak-RAM optimization only, not a correctness fix
    """
    sbwt difference a.sbwt b.sbwt -o ${diff_index} -t ${task.cpus} -v ${low_ram_flag}
    """
}