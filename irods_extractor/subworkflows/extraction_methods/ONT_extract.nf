/*
    ONT WORKFLOWS AND FUNCTIONS
*/

include { BATON_GET_COLLECTION } from '../../modules/baton.nf'

def setONTMetadata(collection_path, linked_metadata) {
    def metadata = [:]
    metadata.irods_path = collection_path
    metadata.read_type = "ont"

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
    | map{irods_item ->
        metamap = [:]
        metamap = setONTMetadata(irods_item.collection, irods_item.avus)
        metamap
    }
    | filter{ meta -> meta.irods_path.contains( params.pre_basecalled ? "pod5": "fastq_pass" ) }
    | set{ meta_ch }

    emit:
    meta_ch
}

workflow ONT_EXTRACT {
    take:
    meta_ch

    main:
    BATON_GET_COLLECTION(meta_ch)
    | set { reads_ch }

    emit:
    reads_ch
}