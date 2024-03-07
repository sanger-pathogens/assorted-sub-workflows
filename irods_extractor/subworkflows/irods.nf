include { COLLATE_FASTQ; COMBINE_FASTQ } from '../modules/samtools.nf'
include { BATON                      } from '../modules/baton.nf'
include { JSON_PREP; JSON_PARSE      } from '../modules/jq.nf'
include { RETRIEVE_CRAM              } from '../modules/retrieve.nf'
include { METADATA                   } from '../modules/metadata_save.nf'

def split_metadata(collection_path, data_obj_name, linked_metadata) {
    metadata = [:]
    metadata.ID = data_obj_name.split("\\.")[0]
    metadata.irods_path = "${collection_path}/${data_obj_name}"
    linked_metadata.each { item ->
        metadata[item.attribute] = item.value
    }
    if (metadata.alt_process){
        metadata.ID = "${metadata.ID}_${metadata.alt_process}"
    }
    return metadata
}

def map_from_multiple(listOfMaps){
    def originMap = listOfMaps.find { it.target == '1' } //select the meta with target == 1 as it is the most complete normally
    def resultMap = [:]
    originMap.each { key, value ->
        if (key == "ID") {
            resultMap[key] = "${value}_total"
        } else if (key == "alignment") {
            resultMap[key] = "Combined"
        } else{
            resultMap[key] = value
        }
    }
    return resultMap
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
        .map{irods_item ->
            metamap = [:]
            metamap = split_metadata(irods_item.collection, irods_item.data_object, irods_item.avus)
            cram_path = metamap.irods_path
            [metamap, cram_path]
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

workflow CRAM_EXTRACT_MERGE {
    
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

    /*
    now we differ from the normal cram extract first we filter out a clean and easy 1 aligned read pair - the golden child
    */
    
    COLLATE_FASTQ.out.fastq_channel
    | filter{ metadata, read_1, read_2 -> (metadata.alignment == '1') }
    | set{ single_channel }

    /*
    We now split out an identifier from the IRODS metadata in the map and append it as a key to use to group
    */

    COLLATE_FASTQ.out.fastq_channel.map{ metaMap, read_1, read_2 ->
        def joinID = "${metaMap.id_run}_${metaMap.lane}#${metaMap.tag_index}"
        tuple(joinID, metaMap, read_1, read_2)
    }.set{ identifer_fastq_ch }

    identifer_fastq_ch.groupTuple().map{ identifier, metadata, read_1, read_2 ->
        tuple(map_from_multiple(metadata), read_1, read_2) //function that makes an amalgam metamap
    }.set{ cleaned_total_reads }

    COMBINE_FASTQ(cleaned_total_reads)

    if (params.cleanup_intermediate_files_irods_extractor) {
        COLLATE_FASTQ.out.remove_channel.flatten()
                .filter(Path)
                .map { it.delete() }
    }

    emit:
    reads_ch = single_channel.mix(COMBINE_FASTQ.out.fastq_channel) // tuple val(meta), path(forward_fastq), path(reverse_fastq
}

workflow IRODS_EXTRACTOR {

    take:
    input_irods_ch // map

    main:

    IRODS_QUERY(input_irods_ch)

    if (params.combine_same_id_crams) {
        /*
        choose your flavor of crams as there are a nice selection
        */
        IRODS_QUERY.out.meta_cram_ch.filter{meta, cram_path -> (meta.ID.endsWith("_human") && meta.alt_process == params.dehumanising_method) || meta.ID.endsWith("_phix") || meta.alignment == '1'}
        .set{ filtered_crams }

        reads_ch = CRAM_EXTRACT_MERGE(filtered_crams)
    } else {
        reads_ch = CRAM_EXTRACT(IRODS_QUERY.out.meta_cram_ch)
    }

    emit:
    reads_ch // tuple val(meta), path(forward_fastq), path(reverse_fastq)
}