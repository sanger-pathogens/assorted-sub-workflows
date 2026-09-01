process THEMISTO2_BUILD {
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

    container "quay.io/sangerpathogens/themisto2:0.0.1"

    // meta.stage disambiguates the candidate rebuild from that same lineage's own
    // THEMISTO2_BUILD_GROUP output -- see ggcat.nf's publishDir for why meta.ID alone can't.
    publishDir mode: 'copy', path: "${params.outdir}/themisto2/${meta.stage ? "${meta.stage}/" : ''}${meta.ID}/build/" 

    input:
    tuple val(meta), path(file_colors_input), path(sbwt_index), path(lcs_index)

    output:
    tuple val(meta), path(index_thm2), emit: index

    script:
    index_thm2 = "index.thm2"
    index_build_params = "--file-colors ${file_colors_input} -o ${index_thm2} -s ${sbwt_index} -l ${lcs_index} -k ${params.color_index_kmer_size} -t ${task.cpus}"

    // User-provided temp storage if given otherwise use tmp workdir (since scratch is enabled).
    // Per-meta.ID subpath -- see sbwt.nf's temp_dir for why (a shared params.temp_dir
    // path would otherwise collide across concurrent species/lineage/candidate tasks).
    if (params.temp_dir) {
        temp_storage_location = "${params.temp_dir}/themisto2/${meta.ID}"
        index_build_params += " --temp-dir ${temp_storage_location}"
    } else {
        temp_storage_location = "\$PWD"
        index_build_params += " --temp-dir ${temp_storage_location}"
    }

    // No memory flag exists in themisto2 build's CLI (that's a v1-only flag) --
    // task.memory is enforced by the scheduler via the mem_* label only.

    """
    mkdir -p ${temp_storage_location}
    sed -i '/^\s*\$/d' "${file_colors_input}"    # Remove blank lines
    themisto2 build ${index_build_params}
    """
}

process THEMISTO2_STATS {
    tag "${meta.ID}"
    label 'cpu_4'
    label 'mem_2' // sized from real test-run peak RSS (<35 MB) -- much lighter than the build, as expected
    label 'time_30m'

    container "quay.io/sangerpathogens/themisto2:0.0.1"

    input:
    tuple val(meta), path(index_thm2)

    output:
    tuple val(meta), path(index_thm2), emit: index

    script:
    // sanity check the index loads and reports counts -- cross-check k-mer/colour
    // counts against Step 04's SBWT index and Step 02's colour-mapping output by eye,
    // this only catches a hard failure to load, not a silent count mismatch.
    """
    themisto2 stats -i ${index_thm2} -t ${task.cpus}
    """
}

process THEMISTO2_EXPORT {
    tag "${meta.ID}"
    label 'cpu_4'
    label 'mem_4'
    label 'time_queue_from_normal'

    container "quay.io/sangerpathogens/themisto2:0.0.1"

    publishDir mode: 'copy', path: "${params.outdir}/themisto2/${meta.stage ? "${meta.stage}/" : ''}${meta.ID}/export/"

    input:
    tuple val(meta), path(index_thm2)

    output:
    tuple val(meta), path("export.unitigs.fa"),    emit: unitigs
    tuple val(meta), path("export.color_sets.txt"), emit: color_sets
    tuple val(meta), path("export.metadata.txt"),  emit: metadata

    script:
    // Export is left uncompressed: single-threaded gzip of color_sets.txt was
    // ~2.3h projected at s_pneumoniae scale, disk has never been the constraint,
    // and the downstream scripts read it directly. (core_catchall_filter.py /
    // lineage_specificity_score.py still transparently accept a .gz if a user
    // supplies a pre-compressed export from elsewhere.)
    //
    // -o is a filename prefix, not a directory -- themisto2 appends its own
    // dot-leading suffixes (.unitigs.fa, .color_sets.txt, .metadata.txt) directly
    // onto it. "-o ." previously produced "..unitigs.fa" etc (prefix "." + suffix
    // ".unitigs.fa"), not the "export.*" names declared in `output:` above.
    """
    themisto2 export -i ${index_thm2} -o export -t ${task.cpus}
    """
}
