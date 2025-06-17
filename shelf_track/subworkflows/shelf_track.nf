include { SHELF_GET_RUN_UUID; SHELF_GET_METHOD_UUID; process SHELF_CREATE_FILE } from '../modules/shelf.nf'

//
// SUBWORKFLOW: Read in study, run, etc. parameters and pull data from iRODS
//
workflow SHELF_TRACK {

    take:
    generic_output // tupple val(meta), path(results)

    main:

    SHELF_GET_RUN_UUID(generic_output)

    SHELF_GET_METHOD_UUID()

    generic_output.join(SHELF_GET_RUN_UUID.out.run_uuid, SHELF_GET_METHOD_UUID.out.method_uuid)
    .set{ file_payload }

    SHELF_CREATE_FILE(file_payload)

    emmit:
    SHELF_CREATE_FILE.out.file_uuid
}