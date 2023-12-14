process BATON {
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'
    maxForks = 10
    container "/software/pathogen/images/baton.simg"
    input:
    path(json_file)

    output:
    path(lane_file), emit: path_channel

    script:
    lane_file="info.json"
    """
    baton-do --file ${json_file} --zone seq > ${lane_file}
    """
}