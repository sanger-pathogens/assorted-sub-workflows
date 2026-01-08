/*
 * SKIP_DOWNLOADED
 *
 * Reusable subworkflow that filters inputs based on presence of
 * existing outputs in an output directory.
 *
 */

workflow FILTER_EXISTING_OUTPUTS {
    take:
    meta_dataobj_ch  // tuple(meta, dataobj); meta is a Map with ID fields, dataobj is the input file/path to be processed 

    main:
    if ("${params.save_method}" == "nested"){
        Channel.fromPath("${params.outdir}/*/${params.preexisting_output_tag}/*${params.existing_output_id_suffix}${params.existing_output_extension}")
        .set{ preexisting_output_path_ch }
    }else{
        Channel.fromPath("${params.outdir}/${params.preexisting_output_tag}/*${params.existing_output_id_suffix}${params.existing_output_extension}")
        .set{ preexisting_output_path_ch }
    }
    preexisting_output_path_ch.toList().map{ preexisting_output_path_list -> 
       no_download = preexisting_output_path_list.size()
       log.info "irods_extractor: ${no_download} data items already exist and won't be downloaded."
    }

    preexisting_output_path_ch.map{ preexisting_output_path ->
        preexisting_output_path.Name.split("${params.existing_output_id_suffix}")[0]
    }
    .collect().map{ [it] }
    .ifEmpty("fresh_run")
    .set{ existing_id }

    meta_dataobj_ch.combine( existing_id )
    .filter { metadata, cram_path, existing -> !(metadata.ID.toString() in existing) }
    .map { it[0,1] }
    .set{ do_not_exist }
    
    do_not_exist.toList().map{ do_not_exist_list ->
        new_downloads = do_not_exist_list.size()
        log.info "irods_extractor: ${new_downloads} data items will be downloaded."
    }

    emit:
    do_not_exist
}
