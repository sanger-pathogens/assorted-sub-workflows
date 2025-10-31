include { PIPELINE_GET_METHOD ;
          PIPELINE_EVENTS_OPEN_BATCH ;
          PIPELINE_EVENTS_CLOSE_BATCH ;
          GATHER_RESULTFILE_INFO ;
          PIPELINE_EVENTS_CREATE_FILE } from '../modules/pipeline_events.nf'


workflow PIPELINE_EVENTS_INIT {
    main:
    PIPELINE_GET_METHOD()
    | PIPELINE_EVENTS_OPEN_BATCH

    batch_uuid = PIPELINE_EVENTS_OPEN_BATCH.out.batch_uuid

    GATHER_RESULTFILE_INFO(PIPELINE_EVENTS_OPEN_BATCH.out.batch_manifest_params, "pipeline_info", "batch_manifest", batch_uuid)
    | PIPELINE_EVENTS_CREATE_FILE

    emit:
    batch_uuid
    batch_manifest_info = PIPELINE_EVENTS_CREATE_FILE.out.created_file_info
}

workflow PIPELINE_EVENTS_END {
    take:
    batch_uuid
    batch_manifest_info
    created_file_infos

    main:
    created_file_infos
    .mix(batch_manifest_info)
    .tap { all_created_file_infos }
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

    PIPELINE_EVENTS_CLOSE_BATCH(batch_uuid, created_file_count)

    emit:
    all_created_file_paths
    all_created_file_infos
    created_file_count
    created_file_count_per_type
}