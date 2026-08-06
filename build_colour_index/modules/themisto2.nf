process THEMISTO_BUILD {
    tag "${meta.ID}"
    label 'cpu_32'
    label 'mem_32'
    label 'time_queue_from_small'

    container "gitlab-registry.internal.sanger.ac.uk/sanger-pathogens/docker-images/themisto2:0.0.1"

    publishDir mode: 'copy', path: "${params.outdir}/themisto2/${meta.ID}/build/"

    input:
    tuple val(meta), path(file_colors_input), path(sbwt_index), path(lcs_index)

    output:
    tuple val(meta), path(index_thm2), emit: index

    script:
    index_thm2 = "index.thm2"
    """
    themisto2 build \\
        --file-colors ${file_colors_input} \\
        -o ${index_thm2} \\
        -s ${sbwt_index} \\
        -l ${lcs_index} \\
        -k ${params.kmer_size} \\
        -t ${task.cpus}

    # sanity check the index loads and reports counts -- cross-check k-mer/colour
    # counts against Step 04's SBWT index and Step 02's colour-mapping output by eye,
    # this only catches a hard failure to load, not a silent count mismatch.
    themisto2 stats -i ${index_thm2} -t ${task.cpus}
    """
}

process THEMISTO_EXPORT {
    tag "${meta.ID}"
    label 'cpu_4'
    label 'mem_16'
    label 'time_queue_from_normal'

    container "gitlab-registry.internal.sanger.ac.uk/sanger-pathogens/docker-images/themisto2:0.0.1"

    publishDir mode: 'copy', path: "${params.outdir}/themisto2/${meta.ID}/export/"

    input:
    tuple val(meta), path(index_thm2)

    output:
    tuple val(meta), path("export.unitigs.fa"),     emit: unitigs
    tuple val(meta), path("export.color_sets.txt*"), emit: color_sets
    tuple val(meta), path("export.metadata.txt"),   emit: metadata

    script:
    // gzip is a separate, deliberately-off-by-default step, not a themisto2 export
    // flag -- single-threaded compression was slow enough at s_pneumoniae scale
    // (~2.3h projected) that the project chose to just leave export uncompressed.
    def gzip_cmd = params.gzip_export ? "gzip -f export.color_sets.txt" : ""
    """
    themisto2 export -i ${index_thm2} -o .
    ${gzip_cmd}
    """
}
