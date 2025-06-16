include { SHELF_GET_RUN_UUID; SHELF_GET_METHOD_UUID; SHELF_PUT_FILE } from '../modules/shelf.nf'

//
// SUBWORKFLOW: Read in study, run, etc. parameters and pull data from iRODS
//
workflow SHELF_TRACK {

    take:
    val(meta)
    path(results)

    main:

    SHELF_GET_RUN_UUID(meta)

    SHELF_GET_METHOD_UUID(results)

    SHELF_PUT_FILE(results, SHELF_GET_RUN_UUID.out.run_uuid, SHELF_GET_METHOD_UUID.out.method_uuid)

    emmit:
    SHELF_PUT_FILE.out.file_uuid
}