process THEMISTO2_BUILD {
    tag "${meta.ID}"
    label 'cpu_32'
    label "mem_8"
    label 'time_queue_from_normal'

    // request /tmp only if /tmp is actually the temp dir (assumes TMPDIR unset)
    if (!params.temp_dir || params.temp_dir.startsWith("/tmp")) {
        label 'request_temp'
    }
    scratch true

    container "quay.io/sangerpathogens/themisto2:0.0.1"

    publishDir mode: 'copy', path: "${params.outdir}/themisto2/${meta.stage ? "${meta.stage}_" : ''}${meta.ID}_build/"

    input:
    tuple val(meta), path(file_colors_input), path(sbwt_index), path(lcs_index)

    output:
    tuple val(meta), path(index_thm2), emit: index

    script:
    index_thm2 = "index.thm2"
    index_build_params = "--file-colors ${file_colors_input} -o ${index_thm2} -s ${sbwt_index} -l ${lcs_index} -k ${params.color_index_kmer_size} -t ${task.cpus}"

    if (params.temp_dir) {
        temp_storage_location = "${params.temp_dir}/themisto2/${meta.ID}"
        index_build_params += " --temp-dir ${temp_storage_location}"
    } else {
        temp_storage_location = "\$PWD"
        index_build_params += " --temp-dir ${temp_storage_location}"
    }


    """
    mkdir -p ${temp_storage_location}
    sed -i '/^\s*\$/d' "${file_colors_input}"    # Remove blank lines
    themisto2 build ${index_build_params}
    """
}

process THEMISTO2_STATS {
    tag "${meta.ID}"
    label 'cpu_4'
    label 'mem_2' // real peak RSS <35 MB
    label 'time_30m'

    container "quay.io/sangerpathogens/themisto2:0.0.1"

    input:
    tuple val(meta), path(index_thm2)

    output:
    tuple val(meta), path(index_thm2), emit: index

    script:
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

    publishDir mode: 'copy', path: "${params.outdir}/themisto2/${meta.stage ? "${meta.stage}_" : ''}${meta.ID}_export/"

    input:
    tuple val(meta), path(index_thm2)

    output:
    tuple val(meta), path("export.unitigs.fa"),    emit: unitigs
    tuple val(meta), path("export.color_sets.txt"), emit: color_sets
    tuple val(meta), path("export.metadata.txt"),  emit: metadata

    script:
    """
    themisto2 export -i ${index_thm2} -o export -t ${task.cpus}
    """
}
