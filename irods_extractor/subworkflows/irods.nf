include { COLLATE_FASTQ; COMBINE_FASTQ  } from '../modules/samtools.nf'
include { BATON                         } from '../modules/baton.nf'
include { JSON_PREP; JSON_PARSE         } from '../modules/jq.nf'
include { RETRIEVE_CRAM                 } from '../modules/retrieve.nf'
include { METADATA as METADATA_QUERIED  } from '../modules/metadata_save.nf'
include { METADATA as METADATA_COMBINED } from '../modules/metadata_save.nf'

def set_metadata(collection_path, data_obj_name, linked_metadata, combine_same_id_crams) {
    def metadata = [:]
    metadata.irods_path = "${collection_path}/${data_obj_name}"
    linked_metadata.each { item ->
        metadata[item.attribute.replaceAll("\\n|\\r", " ")] = item.value
    }

    metadata.subset = "target"
    // target subset remains implicit in ID and file names

    metadata.ID = "${metadata.id_run}_${metadata.lane}${params.lane_plex_sep}${metadata.tag_index}"
    // need to join on 'alt_process' as well, otherwise will combine reads from n different alternative processing options = n x the raw read set
    metadata.ID = !metadata.alt_process ? "${metadata.ID}" : "${metadata.ID}_${metadata.alt_process}"
    if ( combine_same_id_crams ) {
        def slurper = new groovy.json.JsonSlurper()
        def component = slurper.parseText(metadata["component"])
        if ( component.subset ){
            metadata.ID = "${metadata.ID}_${component.subset}"
            metadata.subset = component.subset
        }
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
            metamap = set_metadata(irods_item.collection, irods_item.data_object, irods_item.avus, params.combine_same_id_crams)
            cram_path = metamap.irods_path
            [metamap, cram_path]  }
        .ifEmpty { error("""
            Error: IRODS search returned no data!
            - Please ensure you are logged on with `iinit` and re-run the
            pipeline without `-resume`.
            - If you are logged in, check the process IRODS_QUERY:BATON work
            directory for permission errors and ensure the study exists.
            """)
        }
        .filter{ it[0]["subset"] != "${params.irods_subset_to_skip}" }
        .set{ meta_cram_ch }
        

        if (params.save_metadata) {
            meta_cram_ch.map{metadata_map, path -> metadata_map}
            | collectFile() { map -> [ "lane_metadata.txt", map.collect{it}.join(', ') + '\n' ] }
            | set{ metadata_only }

            metadata_tag = channel.value("irods_queried")
            METADATA_QUERIED(metadata_only, metadata_tag)
        }

        emit:
        meta_cram_ch
}

workflow CRAM_EXTRACT {

    take:
    meta_cram_ch //tuple meta, cram_path

    main:

    if ("${params.save_method}" == "nested"){
        Channel.fromPath("${params.outdir}/*/${params.preexisting_fastq_tag}/*_1.fastq.gz")
        .set{ preexisting_fastq_path_ch }
    }else{
        Channel.fromPath("${params.outdir}/${params.preexisting_fastq_tag}/*_1.fastq.gz")
        .set{ preexisting_fastq_path_ch }
    }
    preexisting_fastq_path_ch.toList().map{ preexisting_fastq_path_list -> 
       no_download = preexisting_fastq_path_list.size()
       log.info "irods_extractor: ${no_download} data items already exist and won't be downloaded."
    }

    preexisting_fastq_path_ch.map{ preexisting_fastq_path ->
        preexisting_fastq_path.Name.split("${params.split_sep_for_ID_from_fastq}")[0]
    }
    .collect().map{ [it] }
    .ifEmpty("fresh_run")
    .set{ existing_id }

    meta_cram_ch.combine( existing_id )
    .filter { metadata, cram_path, existing -> !(metadata.ID.toString() in existing) }
    .map { it[0,1] }
    .set{ do_not_exist }

    do_not_exist.toList().map{ do_not_exist_list -> 
       new_downloads = do_not_exist_list.size()
       log.info "irods_extractor: ${new_downloads} data items will be downloaded."
    }

    RETRIEVE_CRAM(do_not_exist)
    | COLLATE_FASTQ

    if (params.combine_same_id_crams) {
        COLLATE_FASTQ.out.fastq_channel.map{ metaMap, read_1, read_2 ->
            commonid = "${metaMap.id_run}_${metaMap.lane}${params.lane_plex_sep}${metaMap.tag_index}"
            commonid = !metaMap.alt_process ? "${commonid}" : "${commonid}_${metaMap.alt_process}"
            tuple(commonid, metaMap, read_1, read_2)
        }.groupTuple()
        .filter{ it[1].size() >= 2 } // only combine datasets into 'total' subset if more than one subset to start with
        .map{ common_id, metadata_list, read_1_list, read_2_list ->
            tuple(meta_map_for_total_reads(metadata_list), read_1_list.join(' '), read_2_list.join(' ')) // amalgam metamap + concatenated path of read files
        }.filter { it[0] != "none" }
        .set{ gathered_total_reads }

        COMBINE_FASTQ(gathered_total_reads)
        
        COMBINE_FASTQ.out.fastq_channel.mix(COLLATE_FASTQ.out.fastq_channel).set{ reads_ch }

        if (params.save_metadata) {
            COMBINE_FASTQ.out.fastq_channel.map{mdata_map, read1_path, read2_path -> mdata_map}
            | collectFile() { map -> [ "lane_metadata.txt", map.collect{it}.join(', ') + '\n' ] }
            | set{ mdata_only }

            mdata_tag = channel.value("combined_readsets")
            METADATA_COMBINED(mdata_only, mdata_tag)
        }
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