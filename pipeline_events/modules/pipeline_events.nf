

process PIPELINE_GET_METHOD {
    // prpepare Shelf tracking of output files - to apply once for the whole pipeline run
    label 'cpu_1'
    label 'mem_1'
    label 'time_from_queue_small'

    // no input as really we only need to query this once per run based on pipeline own info
    //input:
    //val(pipeline_manifest) // map

    output:
    env(methodurl),  emit: method_url
    env(methodname), emit: method_name // might not match what we've decided to put in Shelf, but good way to provide consistent value or auto-populate
    env(methodshort), emit: method_short // same as above, might be used as a convenient label for searches

    script:
    // relying on manifest scope from main config file 
    // but version field is only populated with real value when running code deployed on farm
    // when running piepline code straight from repo, version field still has template value so taking hard coded value instead
    pipelineurl = workflow.manifest.homePage
    methodshort = (pipelineurl as Path).getSimpleName()
    methodurl = workflow.manifest.version == "{{${methodshort}_version}}" ? "${pipelineurl}" : "${pipelineurl}/-/tree/${workflow.manifest.version}"
    methodname = workflow.manifest.name
    """
    echo "method url: $methodurl; method name: $methodname, method short name: $methodshort"
    """
}


process PIPELINE_EVENTS_OPEN_BATCH {
    label 'cpu_1'
    label 'mem_1'
    label 'time_from_queue_small'

    container 'gitlab-registry.internal.sanger.ac.uk/sanger-pathogens/pipeline-event-api/pipeline-event-api:v1.0.0'

    input:
    val(method_url)
    val(methodname)
    val(methodshort)
    output:
    env(batchuuid),  emit: batch_uuid

    script:
    """
    batchuuid=\$(uuidgen -r)
    send_pipeline_event open --batch_id \${batchuuid} --pipeline_name ${methodname} --pipeline_url ${method_url}
    """
}

process PIPELINE_EVENTS_CLOSE_BATCH {
    label 'cpu_1'
    label 'mem_1'
    label 'time_from_queue_small'

    container 'gitlab-registry.internal.sanger.ac.uk/sanger-pathogens/pipeline-event-api/pipeline-event-api:v1.0.0'

    // no input as really we only need to query this once per run based on pipeline own info
    //input:
    //val(pipeline_manifest) // map

    input:
    val(batchuuid)
    val(filecreatedcount)

    script:
    """
    send_pipeline_event close --batch_id ${batchuuid} --status completed --files_created ${filecreatedcount}
    """
}

process PIPELINE_EVENTS_CREATE_FILE {
    label 'cpu_1'
    label 'mem_1'
    label 'time_from_queue_small'

    container 'gitlab-registry.internal.sanger.ac.uk/sanger-pathogens/pipeline-event-api/pipeline-event-api:v1.0.0'

    input:
    tuple val(meta), path(resultfile)
    val(file_type)
    val(batchuuid)
    val(output_folder)

    output:
    // it would be nice to parse the blob to get the file uuid and at least print it out, maybe all at the end of the pipeline run collecting all blob outputs and doing a bulk print of recorded files
    tuple val(meta.ID), path(resultfile),  emit: created_file_id_path

    script:
    filepath = "${output_folder}/${results}"
    file_outblob = 'shelf_create_file_out.json'
    runid = meta.ID
    """
    filemd5=\$(md5sum ${resultfile})
    send_pipeline_event file --batch_id ${batchuuid} --path ${resultfile} --file_type ${file_type} \\
                                --md5sum ${filemd5} --association RUN --association_id ${runid} \\
                                --username \$USER --group \$(id -gn)
    """

    
}