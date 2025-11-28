include { BEEFEATER                     } from '../modules/beefeater.nf'
include { COLLATE_FASTQ                 } from '../modules/samtools.nf'
include { METADATA as METADATA_QUERIED  } from '../modules/metadata_save.nf'
include { PUBLISH_FASTQ                 } from '../modules/publish.nf'

workflow IRODS_QUERY {
        take:
        input_irods_ch // map

        main:
        BEEFEATER(input_irods_ch)
        | splitJson() //split that sample row into metadata
        | set { meta_file_ch }

        if (params.save_metadata) {
            meta_file_ch
            | collectFile() { map -> [ "lane_metadata.txt", map.collect{it}.join(', ') + '\n' ] }
            | set{ metadata_only }

            METADATA_QUERIED(metadata_only, "irods_queried")
        }

        emit:
        meta_file_ch
}

workflow IRODS_EXTRACTOR {

    take:
    input_irods_ch // map

    main:
    
    meta_file_ch.branch{ meta_map ->
            illumina_to_unpack: meta_map.method == "ILLUMINA"

            ONT: meta_map.method == "ONT"

            other: true
        }
        | set { downloaded_objects }

        CRAM_EXTRACT(downloaded_objects.illumina_to_unpack)

        PUBLISH_FASTQ(downloaded_objects.ONT)

    emit:
    reads_ch = CRAM_EXTRACT.out.reads_ch // tuple val(meta), path(forward_fastq), path(reverse_fastq)
}

workflow CRAM_EXTRACT {

    take:
    meta_cram_ch //tuple meta

    main:

    /*
    if ("${params.save_method}" == "nested"){
        Channel.fromPath("${params.outdir}/*-/${params.preexisting_fastq_tag}/*_1.fastq.gz")
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
    */

    COLLATE_FASTQ(meta_cram_ch )
    | set { reads_ch }

    if (params.cleanup_intermediate_files_irods_extractor) {
        COLLATE_FASTQ.out.remove_channel.flatten()
                .filter(Path)
                .map { it.delete() }
    }

    emit: reads_ch = COLLATE_FASTQ.out.fastq_channel // tuple val(meta), path(forward_fastq), path(reverse_fastq
}