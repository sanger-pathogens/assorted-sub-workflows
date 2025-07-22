process SHELF_GET_RUN_UUID {
    label 'cpu_1'
    label 'mem_1'
    label 'time_from_queue_small'

    //container 'gitlab.internal.sanger.ac.uk/sanger-pathogens/shelf/cli/container_registry/shelf_cli:v0.10.1'

    input:
    tuple val(meta), path(results) // results are not used here but this allows using the same generic input channel for SHELF_CREATE_FILE
    output:
    env(runuuid),  emit: run_uuid

    script:
    """
    module load shelf
    export runuuid=\$(shelf get run -q run.name='${meta.ID}' -H run_uuid | tail -n1)
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
    env(methuuid),  emit: method_uuid

    script:
    // relying on manifest scope from main config file 
    // but version field is only populated with real value when running code deployed on farm
    // when running piepline code straight from repo, version field still has template value so taking hard coded value instead
    pipeline_version = workflow.manifest.version == '{{irods_extractor_version}}' ? 'v3.5.2' : workflow.manifest.version
    pipeline_homepage = workflow.manifest.homePage
    """
    module load shelf
    export methuuid=\$(shelf get method -q url='${pipeline_homepage}/-/tree/${pipeline_version}' -H method_uuid | tail -n1)
    """
}

process SHELF_CREATE_FILE {
    label 'cpu_1'
    label 'mem_1'
    label 'time_from_queue_small'

    //container 'gitlab.internal.sanger.ac.uk/sanger-pathogens/shelf/cli/container_registry/shelf_cli:v0.10.1-rc1'

    input:
    tuple val(meta), path(results)
    val(file_type)
    val(run_uuid)
    val(method_uuid)
    val(output_folder)

    output:
    // it would be nice to parse the blob to get the file uuid and at least print it out, maybe all at the end of the pipeline run collecting all blob outputs and doing a bulk print of recorded files
    tuple val(meta), path(file_outblob),  emit: created_file_blob
    tuple val(meta), env(fileuuid)

    script:
    filepath = "${output_folder}/${results}"
    file_outblob = 'shelf_create_file_out.json'
    """
    module load shelf_staging
    fileuuid=\$(uuidgen -r)
    shelf_staging create file -k file_uuid,path,run_uuid,method_uuid,file_type -v "\${fileuuid},${filepath},${run_uuid},${method_uuid},${file_type}" > ${file_outblob}
    """

    
}

/*
process PARSE_FILE_BLOB {
    label 'cpu_1'
    label 'mem_1'
    label 'time_30m'

    container 'quay.io/sangerpathogens/jq:1.6'

    input:
    tuple val(file_blob)

    output:
    env(file_uuid)

    script:
    'created_file_uuids.txt'
    """
    export file_uuid=\$(jq '{file_uuid}' "${file_blob}")
    """
}
*/