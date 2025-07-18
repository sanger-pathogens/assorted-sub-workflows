process SHELF_GET_RUN_UUID {
    label 'cpu_1'
    label 'mem_1'
    label 'time_from_queue_small'

    //container 'gitlab.internal.sanger.ac.uk/sanger-pathogens/shelf/cli/container_registry/shelf_cli:v0.10.1'

    input:
    tuple val(meta), path(results)

    output:
    tuple val(meta), val(runuuid_out),  emit: run_uuid

    script:
    runuuid_out = 'run_uuid.txt'
    """
    module load shelf/v0.10.1
    shelf get run -q run.name=${meta.ID} -H run_uuid | tail -n1 > $runuuid_out
    """
}

process SHELF_GET_METHOD_UUID {
    label 'cpu_1'
    label 'mem_1'
    label 'time_from_queue_small'

    //container 'gitlab.internal.sanger.ac.uk/sanger-pathogens/shelf/cli/container_registry/shelf_cli:v0.10.1'

    // no input as really we only need to query this once per run based on pipeline own info
    //input:
    //val(pipeline_manifest) // map

    output:
    val(methuuid_out),  emit: method_uuid

    script:
    // relying on manifest scope from main config file but that might not be exported during task
    pipeline_version = workflow.manifest.version == '{{irods_extractor_version}}' ? 'v3.5.2' : workflow.manifest.version
    pipeline_homepage = workflow.manifest.homePage
    methuuid_out = 'meth_uuid.txt'
    """
    module load shelf/v0.10.1
    shelf get method -q url="${pipeline_homepage}/-/tree/${pipeline_version}" -H method_uuid | tail -n1 > $methuuid_out
    """
}

process SHELF_CREATE_FILE {
    label 'cpu_1'
    label 'mem_1'
    label 'time_from_queue_small'

    //container 'gitlab.internal.sanger.ac.uk/sanger-pathogens/shelf/cli/container_registry/shelf_cli:v0.10.1-rc1'

    input:
    tuple val(meta), path(results_file)
    val(file_type)
    val(run_uuid)
    val(method_uuid)
    val(output_folder)

    output:
    val(fileuuid_out),  emit: file_uuid

    script:
    filepath = "${output_folder}/${results_file}"
    fileuuid_out = 'file_uuid.txt'
    """
    module load shelf_staging/v0.10.1-rc1
    shelf_staging create file -k path,run_uuid,method_uuid,file_type -v ${filepath},${run_uuid},${method_uuid},${file_type} | tail -n1 > $fileuuid_out
    """

    
}