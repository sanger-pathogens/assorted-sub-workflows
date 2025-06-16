process SHELF_GET_RUN_UUID {
    label 'cpu_1'
    label 'mem_1'
    label 'time_from_queue_small'

    container 'quay.io/sangerpathogens/shelf:0.0'

    input:
    val(meta)

    output:
    tuple val(meta), env(runuuid),  emit: run_uuid

    script:
    """
    runuuid=\$(shelf get run -q run.name=${meta.ID} -H run_uuid)
    """
}

process SHELF_GET_METHOD_UUID {
    label 'cpu_1'
    label 'mem_1'
    label 'time_from_queue_small'

    container 'quay.io/sangerpathogens/shelf:0.0'

    input:
    path(results)

    output:
    env(methuuid),  emit: method_uuid

    script:
    """
    cd ${projectDir}
    checksum=\$(cat main.nf nextflow.config subworkflows/*.nf assorted-sub-workflows/*/subworkflows/*.nf assorted-sub-workflows/*/subworkflows/*.config modules/*.nf assorted-sub-workflows/*/modules/*.nf | md5sum)
    methuuid=\$(shelf get method -q run.name=\$checksum -H method_uuid)
    """
}

process SHELF_PUT_FILE {
    label 'cpu_1'
    label 'mem_1'
    label 'time_from_queue_small'

    container 'quay.io/sangerpathogens/shelf:0.0'

    input:
    path(results)
    val(run_uuid)
    val(method_uuid)

    output:
    env(fileuuid),  emit: file_uuid

    script:
    """
    cd ${projectDir}
    fileuuid=\$(shelf put )
    """

    
}