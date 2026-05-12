/*
 * SKIP_DOWNLOADED
 *
 * Reusable subworkflow that filters inputs based on presence of
 * existing outputs in an output directory.
 *
 */

def validateNonEmptyString(String name, def val) {
    if (val == null || val.toString().trim() == "") {
        error "Parameter '${name}' must be set to a non-empty string when only_new_input=true."
    }
}

def validateSkipDownloadedParams() {
    validateNonEmptyString('preexisting_output_tag', params.preexisting_output_tag)
    validateNonEmptyString('existing_output_extension', params.existing_output_extension)
    validateNonEmptyString('existing_output_id_suffix', params.existing_output_id_suffix)
}


workflow FILTER_EXISTING_OUTPUTS {
    take:
    meta_dataobj_ch  // tuple(meta, ...); meta is a Map with ID fields and must be the first element; remaining elements are arbitrary

    main:
    validateSkipDownloadedParams()

    if (params.save_method == "nested") {
        Channel.fromPath("${params.outdir}/*/${params.preexisting_output_tag}/*${params.existing_output_id_suffix}${params.existing_output_extension}")
            .set{ preexisting_output_path_ch }
    } else {
        Channel.fromPath("${params.outdir}/${params.preexisting_output_tag}/*${params.existing_output_id_suffix}${params.existing_output_extension}")
            .set{ preexisting_output_path_ch }
    }

    preexisting_output_path_ch.toList().map { preexisting_output_path_list ->
        def no_download = preexisting_output_path_list.size()
        log.info "irods_extractor: ${no_download} data items already exist and won't be downloaded."
    }

    preexisting_output_path_ch
        .map { preexisting_output_path ->
            preexisting_output_path.name.split("${params.existing_output_id_suffix}")[0]
        }
        .collect()
        .map { [it] }
        .ifEmpty("fresh_run")
        .set { existing_id }

    meta_dataobj_ch
        .combine( existing_id )
        .filter { !(it[0].ID.toString() in it[-1]) }
        .map { it[0..-2] }
        .set{ do_not_exist }

    do_not_exist.toList().map { do_not_exist_list ->
        def new_downloads = do_not_exist_list.size()
        log.info "irods_extractor: ${new_downloads} data items will be downloaded."
    }
    emit:
    do_not_exist
}
