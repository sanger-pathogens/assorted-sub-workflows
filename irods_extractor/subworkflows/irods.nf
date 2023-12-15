include { COLLATE_CRAM; FASTQ_FROM_COLLATED_BAM } from '../modules/irods/samtools.nf'
include { BATON } from '../modules/irods/baton.nf'
include { JSON_PREP; JSON_PARSE } from '../modules/irods/jq.nf'
include { RETRIEVE_CRAM } from '../modules/irods/retrieve.nf'

def split_metadata(collection_name, linked_metadata) {
    metadata = [:]
    metadata.ID = collection_name.split("/")[-1].split(".cram")[0]
    linked_metadata.each { item ->
        metadata[item.attribute] = item.value
        }
    metadata
}

workflow IRODS_QUERY {
        take:
        input_irods_ch //tuple study, runid

        main:
        JSON_PREP(input_irods_ch)
        | BATON
        | JSON_PARSE

        JSON_PARSE.out.json_file.splitJson(path: "result").map{collection ->
            meta = [:]
            meta = split_metadata(collection.data_object, collection.avus)
            [meta.ID, meta]
        }.set{ lane_metadata }

        JSON_PARSE.out.paths.splitText().map{ cram_path ->
        ID = cram_path.split("/")[-1].split(".cram")[0]
        [ ID, cram_path ]
        }.set{ cram_ch }

        cram_ch.join(lane_metadata).map{join_identifier, path, meta ->
            [meta, path]
        }.set{ meta_cram_ch }

        emit:
        meta_cram_ch

}

workflow CRAM_EXTRACT {
    
    take:
    input_irods_ch //tuple study, runid

    main:

    Channel.fromPath("${params.outdir}/*#*/raw_fastq/*_1.fastq.gz").map{ raw_fastq_path ->
        ID = raw_fastq_path.simpleName.split("_1")[0]
    }.ifEmpty("fresh_run").set{ existing_id }

    input_irods_ch.combine( existing_id | collect | map{ [it] })
    | filter { meta, cram_path, existing -> !(meta.ID in existing)}
    | map { it[0,1] }
    | set{ do_not_exist }

    RETRIEVE_CRAM(do_not_exist)
    | COLLATE_CRAM
    | FASTQ_FROM_COLLATED_BAM

    FASTQ_FROM_COLLATED_BAM.out.remove_channel.flatten()
            .filter(Path)
            .map { it.delete() }

    emit:
    reads_ch = FASTQ_FROM_COLLATED_BAM.out.fastq_channel
}

workflow IRODS_EXTRACTOR {

    take:
    input_irods_ch //tuple study, runid

    main:
    IRODS_QUERY(input_irods_ch).set{ meta_cram_ch }

    CRAM_EXTRACT(meta_cram_ch)

    emit:
    reads_ch = CRAM_EXTRACT.out
}