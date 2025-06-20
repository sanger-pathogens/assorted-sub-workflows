process SHELF_GET_METHOD_UUID {
    label 'cpu_1'
    label 'mem_1'
    label 'time_from_queue_small'

    container 'quay.io/sangerpathogens/shelf:0.0'

    // no input as really we only need to query this once per run based on pipeline own info

    output:
    env(methuuid),  emit: method_uuid

    script:
    // relying on manifest scope from main config file but that might not be exported during task
    """
    methuuid=\$(shelf get method -q url="${manifest.homePage}/-/tree/${manifest.version}" -H method_uuid)
    """
}

process SHELF_CREATE_COLLECTION {
    label 'cpu_1'
    label 'mem_1'
    label 'time_from_queue_small'

    container 'quay.io/sangerpathogens/shelf:0.0'

    input:
    val(method_uuid)

    output:
    env(coluuid),  emit: col_uuid
   
    script:
    def colname = params.shelf_collection_name ? "${params.shelf_collection_name}" : "${manifest.name} - ${manifest.version} - ${workflow.start}"
    """
    coluuid=\$(shelf create collection -k name,method_uuid -v $colname,$method_uuid)
    """


}

process SHELF_GET_RUN_UUID {
    label 'cpu_1'
    label 'mem_1'
    label 'time_from_queue_small'

    container 'quay.io/sangerpathogens/shelf:0.0'

    input:
    tuple val(meta), path(results)

    output:
    tuple val(meta), env(runuuid),  emit: run_uuid

    script:
    """
    runuuid=\$(shelf get run -q run.name=${meta.ID} -H run_uuid)
    """
}

process SHELF_CREATE_FILE {
    label 'cpu_1'
    label 'mem_1'
    label 'time_from_queue_small'

    container 'gitlab.internal.sanger.ac.uk/sanger-pathogens/shelf/cli/container_registry/shelf_cli:v0.10.1'

    input:
    tuple val(meta), path(results), val(run_uuid), val(method_uuid)

    output:
    env(fileuuid),  emit: file_uuid

    script:
    """
    fileuuid=\$(shelf create file -k run_uuid,method_uuid -v $run_uuid,$method_uuid)
    """

    
}