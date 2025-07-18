include { SHELF_GET_RUN_UUID; SHELF_GET_METHOD_UUID; SHELF_CREATE_FILE } from '../modules/shelf.nf'

//
// SUBWORKFLOW: Prepare registration of pipeline run and register data items in Shelf
//

workflow SHELF_PREPARE {
    // prpepare Shelf tracking of output files - to apply once for the whole pipeline run
    main:
    SHELF_GET_METHOD_UUID()

    SHELF_GET_METHOD_UUID.out.method_uuid
    .set{ method }

    emit:
    method
}


workflow SHELF_TRACK {

    take:
    generic_output      // tupple val(meta), path(results)
    filetype            // val(file_type)
    method              // val(method_uuid)
    outdir              // val(output_folder)

    main:

    SHELF_GET_RUN_UUID(generic_output)

    SHELF_CREATE_FILE(generic_output, filetype, SHELF_GET_RUN_UUID.out.run_uuid, method, outdir)

    emit:
    SHELF_CREATE_FILE.out.file_uuid
}