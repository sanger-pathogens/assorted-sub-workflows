process BATON {
    label 'cpu_1'
    label 'mem_1'
    label 'time_30m'
    maxForks = 10

    container "ghcr.io/wtsi-npg/ub-16.04-baton-irods-4.2.7:5.0.1"
    
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