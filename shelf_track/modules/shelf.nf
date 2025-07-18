process SHELF_GET_RUN_UUID {
    label 'cpu_1'
    label 'mem_1'
    label 'time_from_queue_small'

    //container 'gitlab.internal.sanger.ac.uk/sanger-pathogens/shelf/cli/container_registry/shelf_cli:v0.10.1'

    input:
    tuple val(meta), path(results)

    output:
    tuple val(meta), env(runuuid),  emit: run_uuid

    script:
    """
    module load shelf/v0.10.1
    runuuid=\$(shelf get run -q run.name=${meta.ID} -H run_uuid | tail -n1)
    """
}

process SHELF_GET_METHOD_UUID {
    label 'cpu_1'
    label 'mem_1'
    label 'time_from_queue_small'

    //container 'gitlab.internal.sanger.ac.uk/sanger-pathogens/shelf/cli/container_registry/shelf_cli:v0.10.1'

    // no input as really we only need to query this once per run based on pipeline own info

    output:
    env(methuuid),  emit: method_uuid

    script:
    // relying on manifest scope from main config file but that might not be exported during task
    """
    module load shelf/v0.10.1
    # methuuid=\$(shelf get method -q url="${manifest.homePage}/-/tree/${manifest.version}" -H method_uuid | tail -n1)
    homepage=\$(grep -A10 "^manifest" ${projectDir}/nextflow.config | grep homePage | python3 -c "import sys; print(sys.stdin.readline().split()[-1].strip('\''))")
    version=\$(grep -A10 "^manifest" ${projectDir}/nextflow.config | grep version | python3 -c "import sys; print(sys.stdin.readline().split()[-1].strip('\''))")
    methuuid=\$(shelf get method -q url="\${homepage}/-/tree/\${version}" -H method_uuid | tail -n1)
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
    env(fileuuid),  emit: file_uuid

    script:
    filepath = "${output_folder}/${results_file}"
    """
    module load shelf_staging/v0.10.1-rc1
    fileuuid=\$(shelf_staging create file -k path,run_uuid,method_uuid,file_type -v ${filepath},${run_uuid},${method_uuid},${file_type} | tail -n1)
    """

    
}