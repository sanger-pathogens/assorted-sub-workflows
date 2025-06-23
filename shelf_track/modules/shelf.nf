process SHELF_GET_RUN_UUID {
    label 'cpu_1'
    label 'mem_1'
    label 'time_from_queue_small'

    container 'gitlab.internal.sanger.ac.uk/sanger-pathogens/shelf/cli/container_registry/shelf_cli:v0.10.1'

    input:
    tuple val(meta), path(results)

    output:
    tuple val(meta), env(runuuid),  emit: run_uuid

    script:
    """
    runuuid=\$(shelf get run -q run.name=${meta.ID} -H run_uuid | tail -n1)
    """
}

process SHELF_GET_METHOD_UUID {
    label 'cpu_1'
    label 'mem_1'
    label 'time_from_queue_small'

    container 'gitlab.internal.sanger.ac.uk/sanger-pathogens/shelf/cli/container_registry/shelf_cli:v0.10.1'

    // no input as really we only need to query this once per run based on pipeline own info

    output:
    env(methuuid),  emit: method_uuid

    script:
    // relying on manifest scope from main config file but that might not be exported during task
    """
    methuuid=\$(shelf get method -q url="${manifest.homePage}/-/tree/${manifest.version}" -H method_uuid | tail -n1)
    """
}

process SHELF_CREATE_FILE {
    label 'cpu_1'
    label 'mem_1'
    label 'time_from_queue_small'

    container 'gitlab.internal.sanger.ac.uk/sanger-pathogens/shelf/cli/container_registry/shelf_cli:v0.10.1-rc1'

    input:
    tuple val(meta), path(results_file)
    val(file_type)
    val(run_uuid)
    val(method_uuid)
    val(output_folder)

    output:
    env(fileuuid),  emit: file_uuid

    script:
    filepath = "${output_folder}/${results_file}"
    """
    fileuuid=\$(shelf_staging create file -k path,run_uuid,method_uuid,file_type -v ${filepath},${run_uuid},${method_uuid},${file_type} | tail -n1)
    """

    
}