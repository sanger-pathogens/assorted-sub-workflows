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

    RETRIEVE_CRAM(do_not_exist.transpose())
    | groupTuple
    | set{ grouped_crams }
    
    grouped_crams.branch{ metaInfo, crams ->
        multiple_crams: (crams.size() > 1)
            [ metaInfo, crams ]

        single_cram: (crams.size() == 1)
            [ metaInfo, crams ]
    }.set{ crams_to_merge }

    COMBINE_CRAMS(crams_to_merge.multiple_crams)

    COMBINE_CRAMS.out.merged_cram_ch.mix(crams_to_merge.single_cram).set{ crams_ready_to_collate }
    
    
    COLLATE_FASTQ(crams_ready_to_collate)

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

workflow IRODS_EXTRACTOR_MERGE_CRAMS{
    take:
    input_irods_ch //tuple studyid, runid, laneid, plexid

    main:

    IRODS_QUERY(input_irods_ch)
    
    IRODS_QUERY.out.meta_cram_ch.map{ metatuple, cram ->
        def clean_id = "${metatuple.id_run}_${metatuple.lane}#${metatuple.tag_index}"
        [ clean_id, metatuple, cram ]
    }.set{ clean_id_cram }
    
    clean_id_cram.groupTuple().multiMap{ clean_id, meta, cramlist ->
        def matchingMeta = meta.find { it.ID == clean_id }

        only_target:
            def shortest_path = cramlist.min{ it.length() } //take the shortest path won't have _human for example
            [ matchingMeta, shortest_path ]

        total:
            def total_meta = generate_metadata_from_map(matchingMeta)
            [ total_meta, cramlist ]
    }
    .set{ cram_branch }

    cram_branch.only_target.mix(cram_branch.total).set{ all_crams }
    
    CRAM_EXTRACT(all_crams)

    emit:
    reads_ch = CRAM_EXTRACT.out // tuple val(meta), path(forward_fastq), path(reverse_fastq)
}