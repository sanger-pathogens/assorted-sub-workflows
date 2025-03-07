/*
    ONT WORKFLOWS AND FUNCTIONS
*/

include { RETRIEVE_FASTQ } from '../../modules/retrieve.nf'

def IRODS_ERROR_MSG = """
    Error: IRODS search returned no data!
    - Please ensure you are logged on with `iinit` and re-run the
    pipeline without `-resume`.
    - If you are logged in, check the process IRODS_QUERY:BATON work
    directory for permission errors and ensure the study exists.
"""

def setONTMetadata(collection_path, linked_metadata) {
    def metadata = [:]
    metadata.irods_path = collection_path

    //ONTRUN-xxx or null
    def matcher = metadata.irods_path =~ /ONTRUN-\d+/
    def run_id = matcher ? matcher[0] : null

    metadata.ont_run = run_id 

    linked_metadata.each { item ->
        metadata[item.attribute.replaceAll("\\n|\\r", " ")] = item.value
    }

    metadata.ID = "${run_id}_${metadata.tag_index}"

    return metadata
}

workflow ONT_PARSE {
    take:
    ont_json

    main:
    ont_json
    | filter{ it.text.contains('"attribute": "ont:tag_identifier"') }
    | splitJson(path: "result")
    | map{irods_item ->
        metamap = [:]
        metamap = setONTMetadata(irods_item.collection, irods_item.avus)
        barcode_path = metamap.irods_path
        [metamap, barcode_path]  }
    | filter{ meta, barcode_path -> meta.irods_path.contains("fastq_pass") }
    | ifEmpty { error(IRODS_ERROR_MSG) }
    | set{ meta_barcode_ch }

    emit:
    meta_barcode_ch
}

workflow ONT_EXTRACT {
    take:
    meta_path_ch //tuple meta, barcode_path

    main:
    RETRIEVE_FASTQ(meta_path_ch)
    | set { reads_ch }

    emit:
    reads_ch
}