include { COLLATE_FASTQ; COMBINE_FASTQ } from '../modules/samtools.nf'
include { BATON                      } from '../modules/baton.nf'
include { JSON_PREP; JSON_PARSE      } from '../modules/jq.nf'
include { RETRIEVE_CRAM              } from '../modules/retrieve.nf'
include { METADATA                   } from '../modules/metadata_save.nf'

def set_metadata(collection_path, data_obj_name, linked_metadata) {
    def metadata = [:]
    metadata.irods_path = "${collection_path}/${data_obj_name}"
    linked_metadata.each { item ->
        metadata[item.attribute.replaceAll("\\n|\\r", " ")] = item.value
    }
    metadata.ID = "${metadata.id_run}_${metadata.lane}${params.lane_plex_sep}${metadata.tag_index}"
    // need to join on 'alt_process' as well, otherwise will combine reads from n different alternative processing options = n x the raw read set
    metadata.ID = !metadata.alt_process ? "${metadata.ID}" : "${metadata.ID}_${metadata.alt_process}"
    def slurper = new groovy.json.JsonSlurper()
    def component = slurper.parseText(metadata["component"])
    if ( component.subset ){
        metadata.ID = "${metadata.ID}_${component.subset}"
        metadata.subset = component.subset
    }else{
        metadata.subset = "target"
        // target subset remains implicit in ID and file names
    }
    return metadata
}

def meta_map_for_total_reads(listOfMaps){
    def originMap = listOfMaps.find { it.subset  == "target" } //select the meta with subset field == 'target' as it is the most complete normally and file name is simple.
    if (!originMap){
        // for when no subset in the group is target, e.g. phix subset will be grouped on its own
        return "none"
    }
    def resultMap = [:]
    originMap.each { key, value ->
        if (key == "ID") {
            resultMap[key] = "${value}_total"
        } else {
            resultMap[key] = value
        }
    }
    resultMap.subset = "total"
    resultMap.total_reads = listOfMaps.collect{ it.total_reads.toInteger() }.sum()
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
            metamap = set_metadata(irods_item.collection, irods_item.data_object, irods_item.avus)
            cram_path = metamap.irods_path
            [metamap, cram_path]  }
        .filter{ it[0]["subset"] != "${params.irods_subset_to_skip}" }
        .set{ meta_cram_ch }
        

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

    if (params.combine_same_id_crams) {
        COLLATE_FASTQ.out.fastq_channel.map{ metaMap, read_1, read_2 ->
            commonid = "${metaMap.id_run}_${metaMap.lane}${params.lane_plex_sep}${metaMap.tag_index}"
            commonid = !metaMap.alt_process ? "${commonid}" : "${commonid}_${metaMap.alt_process}"
            tuple(commonid, metaMap, read_1, read_2)
        }.groupTuple()
        .view()
        .filter{ it[1].size >= 2 } // only combine datasets into 'total' subset if more than one subset to start with
        .map{ common_id, metadata_list, read_1_list, read_2_list ->
            tuple(meta_map_for_total_reads(metadata_list), read_1_list.join(' '), read_2_list.join(' ')) // amalgam metamap + concatenated path of read files
        }.filter { it[0] != "none" }
        .set{ gathered_total_reads }

        COMBINE_FASTQ(gathered_total_reads)
        
        COMBINE_FASTQ.out.fastq_channel.mix(COLLATE_FASTQ.out.fastq_channel).set{ reads_ch }
    }else{
        COLLATE_FASTQ.out.fastq_channel.set{ reads_ch }
    }

    if (params.cleanup_intermediate_files_irods_extractor) {
        COLLATE_FASTQ.out.remove_channel.flatten()
                .filter(Path)
                .map { it.delete() }
    }

    emit: reads_ch // tuple val(meta), path(forward_fastq), path(reverse_fastq
}

workflow IRODS_EXTRACTOR {

    take:
    input_irods_ch // map

    main:

    IRODS_QUERY(input_irods_ch)
    | CRAM_EXTRACT

    emit:
    reads_ch = CRAM_EXTRACT.out.reads_ch // tuple val(meta), path(forward_fastq), path(reverse_fastq)
}