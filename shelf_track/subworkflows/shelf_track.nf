include { SHELF_CREATE_COLLECTION; 
          SHELF_GET_METHOD_UUID; 
          SHELF_GET_RUN_UUID; 
          SHELF_CREATE_FILE 
        } from '../modules/shelf.nf'

//
// SUBWORKFLOW: Prepare registration of pipeline run and register data items in Shelf
//
workflow SHELF_PREPARE {
    // prpepare Shelf tracking of output files - to apply once for the whole pipeline run
    main:
    if (params.shelf_track_folder){
        shelf_track_folder_target = Channel.from(params.shelf_track_folder)
    } else {
        shelf_track_folder_target = Channel.from(params.outdir)
    }

    SHELF_GET_METHOD_UUID()

    shelf_track_folder_target.join(SHELF_GET_METHOD_UUID.out.method_uuid)
    | SHELF_CREATE_COLLECTION()

    emmit:
    SHELF_CREATE_COLLECTION.out.col_uuid

}


workflow SHELF_TRACK {
    // enact Shelf tracking of output files - to apply to each data item from the output channel(s) thoughput the pipeline run

    take:
    generic_output // tupple( val(meta), tupple( path(results1) [, path(results2), ...] ))
    collection     // val(col_uuid)

    main:

    SHELF_GET_RUN_UUID(generic_output)

    generic_output.join(SHELF_GET_RUN_UUID.out.run_uuid, collection)
    .set{ file_payload }

    SHELF_CREATE_FILE(file_payload)

    emmit:
    SHELF_CREATE_FILE.out.file_uuid
}

