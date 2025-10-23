import groovy.json.JsonBuilder

process PIPELINE_GET_METHOD {
    // prpepare Shelf tracking of output files - to apply once for the whole pipeline run
    label 'cpu_1'
    label 'mem_1'
    label 'time_from_queue_small'

    // no input as really we only need to query this once per run based on pipeline own info
    //input:
    //val(pipeline_manifest) // map

    output:
    val(methodurl),  emit: method_url
    val(methodname), emit: method_name // might not match what we've decided to put in Shelf, but good way to provide consistent value or auto-populate
    val(methodshort), emit: method_short // same as above, might be used as a convenient label for searches
    val(pipeline_mani_params), emit: pipeline_manifest_params

    script:
    // relying on manifest scope from main config file 
    // but version field is only populated with real value when running code deployed on farm
    // when running piepline code straight from repo, version field still has template value so taking hard coded value instead
    pipelineurl = workflow.manifest.homePage
    methodshort = (pipelineurl as Path).getSimpleName()
    methodurl = workflow.manifest.version == "{{${methodshort}_version}}" ? "${pipelineurl}" : "${pipelineurl}/-/tree/${workflow.manifest.version}"
    methodname = workflow.manifest.name
    pipeline_mani_params = workflow.manifest + params
    """
    echo "method url: $methodurl; method name: $methodname, method short name: $methodshort"
    echo -e "pipeline manifest and run parameters:\n${pipeline_mani_params_json}"
    """
}


process PIPELINE_EVENTS_OPEN_BATCH {
    label 'cpu_1'
    label 'mem_1'
    label 'time_from_queue_small'

    container 'gitlab-registry.internal.sanger.ac.uk/sanger-pathogens/pipeline-event-api/pipeline-event-api:v1.0.2'

    publishDir "${params.outdir}/pipeline_info/", pattern: "*.json", mode: 'copy', overwrite: true

    input:
    val(method_url)
    val(methodname)
    val(methodshort)
    val(pipeline_mani_params)

    output:
    val(batchuuid),  emit: batch_uuid
    path(batch_mani_params_out), emit: batch_manifest_params

    script:
    batchuuid = UUID.randomUUID().toString()
    pipeline_mani_params["batchuuid"] = batchuuid
    pipeline_mani_params_json = new JsonBuilder(pipeline_mani_params).toPrettyString()
    batch_mani_params_out = "pipeline_manifest_run_params_batch_${batchuuid}.json"
    """
    echo "${pipeline_mani_params_json}" > ${pipeline_mani_params_out}
    send_pipeline_event open --batch_id ${batchuuid} --pipeline_name ${methodname} --pipeline_url ${method_url}
    """
}

process PIPELINE_EVENTS_CLOSE_BATCH {
    label 'cpu_1'
    label 'mem_1'
    label 'time_from_queue_small'

    container 'gitlab-registry.internal.sanger.ac.uk/sanger-pathogens/pipeline-event-api/pipeline-event-api:v1.0.2'

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

    container 'gitlab-registry.internal.sanger.ac.uk/sanger-pathogens/pipeline-event-api/pipeline-event-api:v1.0.2'

    input:
    tuple val(meta), path(resultfileWorkPath), val(resultfilePublishedDir) // val(resultfilePublishedDir), not path() so no to stage folder
    val(file_type)
    val(batchuuid)

    output:
    tuple val(meta.ID), val(resultfilePublishedFullPath),  emit: created_file_id_path // val(resultfilePublishedFullPath), not path() so no to stage the file that's outside the work folder
    val(file_type)

    script:
    runid = meta.ID
    resultfileName = resultfileWorkPath.name.toString()
    resultfilePublishedFullPath = "${resultfilePublishedDir}/${resultfileName}"
    """
    filemd5=\$(md5sum ${resultfileWorkPath} | cut -d' ' -f1)
    send_pipeline_event file --batch_id ${batchuuid} --path ${resultfilePublishedFullPath} --file_type ${file_type} \\
                                --md5sum \${filemd5} --association RUN --association_id ${runid} \\
                                --username \$(id -un) --group \$(id -gn)
    """

    
}

process PIPELINE_EVENTS_CREATE_BATCH_MANIFEST_FILE {
    label 'cpu_1'
    label 'mem_1'
    label 'time_from_queue_small'

    container 'gitlab-registry.internal.sanger.ac.uk/sanger-pathogens/pipeline-event-api/pipeline-event-api:v1.0.2'

    input:
    val(batchuuid)
    path(batch_manifest_params)

    output:
    tuple val(meta.ID), val(resultfilePublishedFullPath),  emit: created_file_path // val(resultfilePublishedFullPath), not path() so no to stage the file that's outside the work folder
    val(file_type)

    script:
    batchManifestfileName = batch_manifest_params.name.toString()
    batchManifestfilePublishedFullPath = "${params.outdir}/pipeline_info/${resultfileName}"
    """
    filemd5=\$(md5sum ${resultfileWorkPath} | cut -d' ' -f1)
    send_pipeline_event file --batch_id ${batchuuid} --path ${resultfilePublishedFullPath} --file_type "batch_manifest" \\
                                --md5sum \${filemd5} \\
                                --username \$(id -un) --group \$(id -gn)
    """

    
}

workflow PIPELINE_EVENTS_INIT {

    main:
    PIPELINE_GET_METHOD() 
    | PIPELINE_EVENTS_OPEN_BATCH
    | PIPELINE_EVENTS_CREATE_BATCH_MANIFEST_FILE 

    emit:
    batch_uuid = PIPELINE_EVENTS_OPEN_BATCH.out.batch_uuid
    batch_manifest = PIPELINE_EVENTS_CREATE_BATCH_MANIFEST_FILE.out.created_file_path
}


workflow PIPELINE_EVENTS_END {

    take:
    batch_uuid
    batch_manifest
    created_file_id_paths

    main:
    created_file_id_paths
    .map{ runid, filepath -> filepath }
    .mix(batch_manifest)
    .count()
    .set{ created_file_count }

    PIPELINE_EVENTS_CLOSE_BATCH(batch_uuid, created_file_id_paths)

    emit:
    created_file_id_paths
}