import groovy.json.JsonBuilder;
import java.nio.file.LinkOption;
import java.lang.reflect.Field as ReflectField;

process PIPELINE_GET_METHOD {
    // prpepare Shelf tracking of output files - to apply once for the whole pipeline run
    label 'cpu_1'
    label 'mem_1'
    label 'time_queue_from_small'

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
    // workflow.manifest does not have a toMap() method, so building map manually; need to cast to ReflectField to avoid conflict with groovy Field
    Map pipeline_manifest = [:]
    ReflectField[] fields = workflow.manifest.getClass().getFields();
    for(ReflectField f : fields){
        Object v = f.get(workflow.manifest);
        pipeline_manifest[f.getName()] = v;
    }
    Map pipeline_mani_params = [:]
    pipeline_mani_params["pipeline_manifest"] = pipeline_manifest
    pipeline_mani_params["params"] = params as Map
    """
    echo "method url: $methodurl; method name: $methodname, method short name: $methodshort"
    """
}


process PIPELINE_EVENTS_OPEN_BATCH {
    label 'cpu_1'
    label 'mem_1'
    label 'time_queue_from_small'

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
    // Map batch_mani_params = params as Map
    batch_mani_params["batchuuid"] = batchuuid
    batch_mani_params_json = new JsonBuilder(batch_mani_params).toPrettyString()
    batch_mani_params_out = "pipeline_manifest_run_params_batch_${batchuuid}.json"
    """
    echo "${batch_mani_params_json}" > ${batch_mani_params_out}
    send_pipeline_event open --batch_id ${batchuuid} --pipeline_name ${methodname} --pipeline_url ${method_url}
    """
}

process PIPELINE_EVENTS_CLOSE_BATCH {
    label 'cpu_1'
    label 'mem_1'
    label 'time_queue_from_small'

    container "${params.pipeline_events_container}"

    input:
    val(batchuuid)
    val(filecreatedcount)

    script:
    """
    send_pipeline_event close --batch_id ${batchuuid} --status completed --files_created ${filecreatedcount}
    """
}

process GET_RESULTFILE_PATH {
    label 'cpu_1'
    label 'mem_1'
    label 'local'

    input:
    tuple val(meta), path(resultfileWorkPath)
    val(outputfoldertag)

    output:
    tuple val(meta), path(resultfileWorkPath), val(resultfilePublishedDirAbsPath),  emit: resultfile_publisheddir

    script:
    if (("${params.save_method}" == "flat" & "${outputfoldertag}".endsWith("fastqs")) | (meta.ID == null)){
        resultfilePublishedDirRelPath = "${params.outdir}/${outputfoldertag}"
    } else {
        resultfilePublishedDirRelPath = "${params.outdir}/${meta.ID}/${outputfoldertag}"
    }
    if (!params.pipeline_events_follow_links) {
        resultfilePublishedDirAbsPath = file(resultfilePublishedDirRelPath).toRealPath(LinkOption.NOFOLLOW_LINKS).toString()
    } else {
        resultfilePublishedDirAbsPath = file(resultfilePublishedDirRelPath).toAbsolutePath().toString()
    }
    """
    echo "file path to track: ${resultfilePublishedDirAbsPath}"
    """
}

process PIPELINE_EVENTS_CREATE_FILE {
    label 'cpu_1'
    label 'mem_1'
    label 'time_queue_from_small'

    container "${params.pipeline_events_container}"

    input:
    tuple val(meta), path(resultfileWorkPath), val(resultfilePublishedDir) // val(resultfilePublishedDir), not path() so not to stage folder
    val(file_type)
    val(batchuuid)

    output:
    tuple val(outid), val(resultfilePublishedFullPath), val(file_type),  emit: created_file_id_path // val(resultfilePublishedFullPath), not path() so not to stage the file that's outside the work folder

    script:
    runid = meta.ID
    outid = (runid == null) ? batchuuid : runid
    resultfileName = resultfileWorkPath.name.toString()
    resultfilePublishedFullPath = "${resultfilePublishedDir}/${resultfileName}"
    runassociationOptString = (file_type == "batch_manifest") ? "" : "--association RUN --association_id ${runid}"
    """
    filemd5=\$(md5sum ${resultfileWorkPath} | cut -d' ' -f1)
    send_pipeline_event file --batch_id ${batchuuid} --path ${resultfilePublishedFullPath} --file_type ${file_type} \\
                                --md5sum \${filemd5} ${runassociationOptString} \\
                                --username \$(id -un) --group \$(id -gn)
    """

    
}

// process PIPELINE_EVENTS_CREATE_BATCH_MANIFEST_FILE {
//     label 'cpu_1'
//     label 'mem_1'
//     label 'time_queue_from_small'

//     container "${params.pipeline_events_container}"

//     input:
//     tuple val(batch_mani_params), path(batchManifestfileWorkPath), val(batchManifestfilePublishedDir) // val(batchManifestfilePublishedDir), not path() so not to stage folder
//     val(batchuuid)

//     output:
//     tuple val(batchuuid), val(batchManifestfilePublishedFullPath), val(file_type),  emit: created_file_path // val(batchManifestfilePublishedFullPath), not path() so not to stage the file that's outside the work folder

//     script:
//     batchManifestfileName = batchManifestfileWorkPath.name.toString()
//     batchManifestfilePublishedFullPath = "${batchManifestfilePublishedDir}/${batchManifestfileName}"
//     """
//     filemd5=\$(md5sum ${batchManifestfileWorkPath} | cut -d' ' -f1)
//     send_pipeline_event file --batch_id ${batchuuid} --path ${batchManifestfilePublishedFullPath} --file_type "batch_manifest" \\
//                                 --md5sum \${filemd5} \\
//                                 --username \$(id -un) --group \$(id -gn)
//     """
// }

workflow PIPELINE_EVENTS_INIT {

    main:
    PIPELINE_GET_METHOD() 
    | PIPELINE_EVENTS_OPEN_BATCH
    
    GET_RESULTFILE_PATH(PIPELINE_EVENTS_OPEN_BATCH.out.batch_manifest_params, "pipeline_info")

    PIPELINE_EVENTS_CREATE_FILE(GET_RESULTFILE_PATH.out.resultfile_publisheddir, "batch_manifest", PIPELINE_EVENTS_OPEN_BATCH.out.batch_uuid)

    emit:
    batch_uuid = PIPELINE_EVENTS_OPEN_BATCH.out.batch_uuid
    batch_manifest = PIPELINE_EVENTS_CREATE_FILE.out.created_file_path
}


workflow PIPELINE_EVENTS_END {

    take:
    batch_uuid
    batch_manifest
    created_file_id_paths

    main:
    created_file_id_paths
    .map{ runid, filepath, filetype -> filepath }
    .mix(batch_manifest)
    .set { all_created_file_paths }

    all_created_file_paths
    .count()
    .set{ created_file_count }

    PIPELINE_EVENTS_CLOSE_BATCH(batch_uuid, created_file_count)

    emit:
    all_created_file_paths
}