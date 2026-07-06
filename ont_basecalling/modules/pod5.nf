process CONVERT_FAST5_TO_POD5 {
    label 'cpu_8'
    label 'mem_8'
    label 'time_30m'

    container 'quay.io/sangerpathogens/pod5:0.3.6'
    
    input:
    tuple val(meta), path(fast5s)

    output:
    tuple val(meta), path(final_name), emit: pod5_ch

    script:
    final_name = "${meta.ID}_converted.pod5"
    """
    pod5 convert fast5 --output ${final_name} -t ${task.cpus} ${fast5s}
    """
}

process MERGE_POD5 {
    label 'cpu_1'
    label 'mem_1'
    label 'time_30m'

    container 'quay.io/sangerpathogens/pod5:0.3.6'
    
    input:
    tuple val(meta), path(pod5s)

    output:
    tuple val(meta), path(final_name), emit: pod5_ch

    script:
    final_name = "${meta.ID}_merged.pod5"
    """
    pod5 merge --output ${final_name} -t ${task.cpus} ${pod5s}
    """
}

process CHECK_POD5_CHEMISTRY {
    label 'cpu_1'
    label 'mem_1'
    label 'time_30m'

    container 'quay.io/sangerpathogens/pod5:0.3.6'
    
    input:
    tuple val(meta), path(pod5)

    output:
    tuple val(meta), path(output_json), path(pod5), emit: pod5_ch

    script:
    output_json = "${meta.ID}_read_info.json"

    command = "${projectDir}/assorted-sub-workflows/ont_basecalling/bin/check_pod5_chemistry.py"
    """
    ${command} ${pod5} -o ${output_json}
    """
}