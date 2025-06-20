include { SHELF_GET_RUN_UUID; SHELF_GET_METHOD_UUID; process SHELF_CREATE_FILE } from '../modules/shelf.nf'

//
// SUBWORKFLOW: Prepare registration of pipeline run and register data items in Shelf
//

workflow SHELF_PREPARE {
    // prpepare Shelf tracking of output files - to apply once for the whole pipeline run
    main:
    SHELF_GET_METHOD_UUID()

    emmit:
    SHELF_GET_METHOD_UUID.out.method_uuid
}


workflow SHELF_TRACK {

    take:
    generic_output // tupple val(meta), path(results)
    method         // val(method_uuid)

    main:

    SHELF_GET_RUN_UUID(generic_output)

    generic_output.join(SHELF_GET_RUN_UUID.out.run_uuid, method)
    | SHELF_CREATE_FILE

    emmit:
    SHELF_CREATE_FILE.out.file_uuid
}