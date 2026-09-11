process THEMISTO2_BUILD {
    tag "${meta.ID}"
    label 'cpu_32'
    label 'mem_16'
    label 'time_queue_from_normal'

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
    sed -i '/^\s*\$/d' "${file_colors_input}"
    themisto2 build ${index_build_params}
    """
}

process THEMISTO2_STATS {
    tag "${meta.ID}"
    label 'cpu_4'
    label 'time_30m'
    label 'mem_16'

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
    label 'mem_16'
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

process THEMISTO2_ATB_PSEUDOALIGN {
    tag "${meta.ID}"
    label 'cpu_4'
    label 'mem_350'
    label 'time_queue_from_normal'

    container "quay.io/sangerpathogens/themisto2:0.0.1"

    publishDir mode: 'copy', path: "${params.outdir}/atb_cross_species/${meta.ID}/", enabled: params.publish_intermediate

    input:
    tuple val(meta), path(candidate_fasta)
    path(atb_index)

    output:
    tuple val(meta), path(jsonl), emit: jsonl

    script:
    jsonl = "${meta.ID}_atb_pseudoalign.jsonl"
    """
    echo "${candidate_fasta}" > querylist.txt
    echo "${jsonl}"           > outlist.txt

    themisto2 threshold-pseudoalign -i ${atb_index} \\
        --query-list querylist.txt \\
        --query-output-list outlist.txt \\
        --threshold 0.01 --denominator all --min-hits 1 --report-hit-counts \\
        -t ${task.cpus}
    """
}
