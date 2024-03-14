include { COLLATE_FASTQ            } from '../modules/samtools.nf'
include { BATON                    } from '../modules/baton.nf'
include { JSON_PREP; JSON_PARSE    } from '../modules/jq.nf'
include { RETRIEVE_CRAM            } from '../modules/retrieve.nf'
include { METADATA                 } from '../modules/metadata_save.nf'

def dataObj_path_to_ID(dataObj_path) {
    split_dataObj_dir = dataObj_path.split("/")
    dataObj_simplename = split_dataObj_dir[-1].split(".")[0]
    // output file renaming logic relies on the iRODS folder structure, which is not perfect, but can't think of anything else there
    ID = { split_dataObj_dir[-2].startsWith('plex') ? dataObj_simplename : "${split_dataObj_dir[-2]}_${dataObj_simplename}" }
    return ID
}

def split_metadata(irodsdatapath, linked_metadata) {
    metadata = [:]
    metadata.ID = dataObj_path_to_ID(irodsdatapath)
    linked_metadata.each { item ->
        metadata[item.attribute] = item.value
        }
    return metadata
}

workflow IRODS_QUERY {
        take:
        input_irods_ch // map

        main:
        JSON_PREP(input_irods_ch)
        | BATON
        | JSON_PARSE

        JSON_PARSE.out.json_file
        .filter{ it.text.contains('"attribute": "alignment"') }
        .splitJson(path: "result")
        .map{collection ->
            metaparse = [:]
            metaparse = split_metadata(collection.data_object, collection.avus)
            [metaparse.ID, metaparse]
        }.set{ lane_metadata }


        JSON_PARSE.out.paths.splitText().map{ cram_path ->
            ID = dataObj_path_to_ID(cram_path)
            [ ID, cram_path ]
        }.set{ cram_ch }

        cram_ch.join(lane_metadata).map{join_identifier, path, metamap ->
            [metamap, path]
        }.set{ meta_cram_ch }

        if (params.save_metadata) {
            meta_cram_ch.map{metadata_map, path -> metadata_map}
            | collectFile() { map -> [ "lane_metadata.txt", map.collect{it}.join(', ') + '\n' ] }
            | set{ metadata_only }

            METADATA(metadata_only)
        }

        emit:
        meta_cram_ch

}

workflow CRAM_EXTRACT {
    
    take:
    meta_cram_ch //tuple meta, cram_path

    main:

    Channel.fromPath("${params.outdir}/*/${params.preexisting_fastq_tag}/*_1.fastq.gz").map{ preexisting_fastq_path ->
        ID = preexisting_fastq_path.Name.split("${params.split_sep_for_ID_from_fastq}")[0]
    }.ifEmpty("fresh_run").set{ existing_id }

    meta_cram_ch.combine( existing_id | collect | map{ [it] })
    | filter { metadata, cram_path, existing -> !(metadata.ID in existing)}
    | map { it[0,1] }
    | set{ do_not_exist }

    RETRIEVE_CRAM(do_not_exist)
    | COLLATE_FASTQ

    if (params.cleanup_intermediate_files_irods_extractor) {
        COLLATE_FASTQ.out.remove_channel.flatten()
                .filter(Path)
                .map { it.delete() }
    }
    emit:
    reads_ch = COLLATE_FASTQ.out.fastq_channel // tuple val(meta), path(forward_fastq), path(reverse_fastq)
}

workflow IRODS_EXTRACTOR {

    take:
    input_irods_ch // map

    main:

    IRODS_QUERY(input_irods_ch)
    | CRAM_EXTRACT

    emit:
    reads_ch = CRAM_EXTRACT.out // tuple val(meta), path(forward_fastq), path(reverse_fastq)
}
