include { COLLATE_FASTQ; COMBINE_CRAMS } from '../modules/samtools.nf'
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

def generate_metadata_from_map(orginMap){
    def resultMap = [:]
    orginMap.each { key, value ->
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
    /*
    Get all the crams first so we can focus on getting the crams packaged as a total unit
    */

    meta_cram_ch.multiMap{ metaTuple, cram ->
        def clean_id = "${metaTuple.id_run}_${metaTuple.lane}#${metaTuple.tag_index}"
        metadata: tuple(clean_id, metaTuple) //no cram for now as it will be grabbed out of the pile later

        crams: tuple(clean_id, cram)
    }.set{ seperated }
    
    Channel.fromPath("${params.outdir}/*/${params.preexisting_fastq_tag}/*_1.fastq.gz").map{ preexisting_fastq_path ->
        ID = preexisting_fastq_path.Name.split("${params.split_sep_for_ID_from_fastq}")[0]
    }.ifEmpty("fresh_run").set{ existing_id }

    seperated.crams.map{ identifer, cram ->
        def metaMap = [:]
        metaMap.ID = identifer
        [ metaMap, cram ]
    }.combine( existing_id | collect | map{ [it] })
    | filter { metadata, cram_path, existing -> !(metadata.ID in existing)}
    | map { it[0,1] }
    | set{ do_not_exist }

    RETRIEVE_CRAM(do_not_exist).map{ metatuple, cram -> tuple(metatuple.ID, cram)}.set { downloaded_crams }

    /*
    now we have crams in the form of [ metatuple.ID, metatuple, all_crams ]
    */

    seperated.metadata.join(downloaded_crams).groupTuple().map{ clean_id, metaMap, cramList ->
            tuple(metaMap, cramList)
    }.set{ ready_to_sort_ch }
    
    ready_to_sort_ch.multiMap{ metaTuples, crams ->
        def shortest_path = crams.min{ it.getFileName().toString().length() } //take the shortest path won't have _human for example
        def matchingMeta = metaTuples.find { it.alignment == '1' } //select the alignment 1 meta
        single: tuple(matchingMeta, shortest_path)

        total: tuple(generate_metadata_from_map(matchingMeta), crams)
    }.set{ sorted_ch }

    COMBINE_CRAMS(sorted_ch.total)

    COMBINE_CRAMS.out.merged_cram_ch.mix(sorted_ch.single).set{ crams_ready_to_collate }
    
    COLLATE_FASTQ(crams_ready_to_collate)

    if (params.cleanup_intermediate_files_irods_extractor) {
        COLLATE_FASTQ.out.remove_channel.flatten()
                .filter(Path)
                .map { it.delete() }
    }

    emit:
    reads_ch = COLLATE_FASTQ.out.fastq_channel // tuple val(meta), path(forward_fastq), path(reverse_fastq
}

workflow IRODS_EXTRACTOR {

    take:
    input_irods_ch // map

    main:

    IRODS_QUERY(input_irods_ch)

    if (params.combine_same_id_crams) {
        reads_ch = CRAM_EXTRACT_MERGE(IRODS_QUERY.out.meta_cram_ch)
    } else {
        reads_ch = CRAM_EXTRACT(IRODS_QUERY.out.meta_cram_ch)
    }

    emit:
    reads_ch // tuple val(meta), path(forward_fastq), path(reverse_fastq)
}