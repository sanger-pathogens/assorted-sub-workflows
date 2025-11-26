include { PIPELINE_GET_METHOD ;
          PIPELINE_EVENTS_OPEN_BATCH ;
          PIPELINE_EVENTS_CLOSE_BATCH ;
          GATHER_RESULTFILE_INFO ;
          PIPELINE_EVENTS_CREATE_FILE } from '../modules/pipeline_events.nf'


workflow PIPELINE_EVENTS_INIT {
    take:
    files_to_track

    main:
    PIPELINE_GET_METHOD()

    files_to_track
    | first()
    | map { metamap, filepath -> true }
    | ifEmpty(false)
    | set { files_created }

    PIPELINE_EVENTS_OPEN_BATCH(
        PIPELINE_GET_METHOD.out.method_url,
        PIPELINE_GET_METHOD.out.method_name,
        PIPELINE_GET_METHOD.out.method_short,
        PIPELINE_GET_METHOD.out.pipeline_manifest_params,
        files_created
    )

    batch_id = PIPELINE_EVENTS_OPEN_BATCH.out.batch_id

    GATHER_RESULTFILE_INFO(PIPELINE_EVENTS_OPEN_BATCH.out.batch_manifest_params, "pipeline_info", "batch_manifest", batch_id)

    if (params.track_batch_metadata) {
        GATHER_RESULTFILE_INFO.out.file_info
        | PIPELINE_EVENTS_CREATE_FILE
        | set { batch_manifest_info }
    } else {
        Channel.empty()
        | set { batch_manifest_info }  // default empty channel if not associating metadata
    }

    emit:
    batch_id
    batch_manifest_info
}

workflow PIPELINE_EVENTS_END {
    take:
    batch_id
    batch_manifest_info
    created_file_infos

    main:
    if (params.associate_batch_metadata) {
        created_file_infos
        .mix(batch_manifest_info)
        .set { all_created_file_infos }
    } else {
        created_file_infos
        .set { all_created_file_infos }
    }

    all_created_file_infos
    .map { runid, filepath, filetype -> filepath }
    .set { all_created_file_paths }

    all_created_file_paths
    .count()
    .set { created_file_count }

    created_file_count
    .subscribe{ file_count -> log.info("Total count of file tracked in Pipeline Events Database: ${file_count}") }

    all_created_file_infos
    .groupTuple(by: 2) // group by file_type
    .map { ids, file_id_paths, file_type -> [file_type, file_id_paths.size()] }
    .set { created_file_count_per_type }

    created_file_count_per_type
    .subscribe { file_count_per_type -> log.info("Count of file tracked in Pipeline Events Database per file type: ${file_count_per_type}") }

    PIPELINE_EVENTS_CLOSE_BATCH(batch_id, created_file_count)

    emit:
    all_created_file_paths
    all_created_file_infos
    created_file_count
    created_file_count_per_type
}