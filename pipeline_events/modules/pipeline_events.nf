import groovy.json.JsonBuilder;

process PIPELINE_GET_METHOD {
    // prpepare Shelf tracking of output files - to apply once for the whole pipeline run
    label 'cpu_1'
    label 'mem_1'
    label 'local'
    cache false

    // no input as really we only need to query this once per run based on pipeline own info

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
    // Collect pipeline metadata into a map
    pipeline_manifest = workflow.manifest.toMap()
    pipeline_mani_params = [:] // https://www.nextflow.io/docs/latest/reference/syntax.html#variable-declaration "Variables declared in the process script, exec, and stub sections exist only in their respective section, with one exception – variables declared without the def keyword also exist in the output section"
    pipeline_mani_params["pipeline_manifest"] = pipeline_manifest
    pipeline_mani_params["params"] = params as LinkedHashMap
    pipeline_mani_params["PWD"] =  System.getenv("PWD")
    """
    echo "method url: $methodurl; method name: $methodname, method short name: $methodshort"
    """
}


process PIPELINE_EVENTS_OPEN_BATCH {
    label 'cpu_1'
    label 'mem_1'
    label 'time_queue_from_small'
    cache false

    container "${params.pipeline_events_container}"

    publishDir "${params.outdir}/pipeline_info/", pattern: "*.json", mode: 'copy', overwrite: true

    input:
    val(method_url)
    val(methodname)
    val(methodshort)
    val(batch_mani_params)

    output:
    tuple val(batch_mani_params), path(batch_mani_params_out), emit: batch_manifest_params
    val(batchuuid),  emit: batch_uuid

    script:
    batchuuid = UUID.randomUUID().toString()
    batch_mani_params["batchuuid"] = batchuuid
    batch_mani_params_json = new JsonBuilder(batch_mani_params).toPrettyString()
    batch_mani_params_out = "pipeline_manifest_run_params_batch_${batchuuid}.json"
    """
    echo '${batch_mani_params_json}' > ${batch_mani_params_out}
    send_pipeline_event open --batch_id ${batchuuid} --pipeline_name ${methodname} --pipeline_url ${method_url}
    """
}

process PIPELINE_EVENTS_CLOSE_BATCH {
    label 'cpu_1'
    label 'mem_1'
    label 'time_queue_from_small'
    cache false

    container "${params.pipeline_events_container}"

    input:
    val(batchuuid)
    val(filecreatedcount)

    script:
    """
    send_pipeline_event close --batch_id ${batchuuid} --status completed --files_created ${filecreatedcount}
    """
}

process GATHER_RESULTFILE_INFO {
    label 'cpu_1'
    label 'mem_1'
    label 'local'
    cache false

    input:
    tuple val(meta), path(resultfileWorkPath)
    val(outputfoldertag)
    val(file_type)
    val(batchuuid)

    output:
    tuple val(meta), path(resultfileWorkPath), val(resultfilePublishedDirAbsPath), val(file_type), val(batchuuid), emit: file_info

    // could be exec block here maybe, given shell script is there purely to avoid returning last declared variable
    script:
    if (!params.pipeline_events_follow_links) {
        pwd = file(System.getenv("PWD"))
        outDirAbsPath = pwd.resolve(params.outdir).normalize().toString()
    } else {
        outDirAbsPath = file(params.outdir).toAbsolutePath().toString()
    }
    if (("${params.save_method}" == "flat" & "${outputfoldertag}".endsWith("fastqs")) | (meta.ID == null)){
        resultfilePublishedDirAbsPath = "${outDirAbsPath}/${outputfoldertag}"
    } else {
        resultfilePublishedDirAbsPath = "${outDirAbsPath}/${meta.ID}/${outputfoldertag}"
    }
    """
    echo "file path to track: ${resultfilePublishedDirAbsPath}"
    """
}

process PIPELINE_EVENTS_CREATE_FILE {
    label 'cpu_1'
    label 'mem_1'
    label 'time_queue_from_small'
    cache false

    container "${params.pipeline_events_container}"

    input:
    tuple val(meta), path(resultfileWorkPath), val(resultfilePublishedDir), val(file_type), val(batchuuid) // val(resultfilePublishedDir), not path() so not to stage folder

    output:
    tuple val(outid), val(resultfilePublishedFullPath), val(file_type),  emit: created_file_info // val(resultfilePublishedFullPath), not path() so not to stage the file that's outside the work folder

    script:
    runid = meta.ID
    outid = (runid == null) ? batchuuid : runid
    resultfileName = resultfileWorkPath.name.toString()
    resultfilePublishedFullPath = "${resultfilePublishedDir}/${resultfileName}"
    // runassociationOptString = (file_type == "batch_manifest") ? "" : "--association RUN --association_id ${runid}"
    runassociationOptString = "--association RUN --association_id ${runid}" // in waiting for the API to allow not using these options for batch_manifest type
    """
    filemd5=\$(md5sum ${resultfileWorkPath} | cut -d' ' -f1)
    send_pipeline_event file --batch_id ${batchuuid} --path ${resultfilePublishedFullPath} --file_type ${file_type} \\
                                --md5sum \${filemd5} ${runassociationOptString} \\
                                --username \$(id -un) --group \$(id -gn)
    """

    
}
