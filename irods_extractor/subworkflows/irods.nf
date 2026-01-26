include { BEEFEATER                     } from '../modules/beefeater.nf'
include { COLLATE_FASTQ                 } from '../modules/samtools.nf'
include { METADATA as METADATA_QUERIED  } from '../modules/metadata_save.nf'
include { PUBLISH_FASTQ                 } from '../modules/publish.nf'

workflow IRODS_QUERY {
        main:
        BEEFEATER()
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
    
    input_irods_ch.branch{ meta_map ->
            illumina_to_unpack: meta_map.Platform == "ILLUMINA"

            ONT: meta_map.Platform == "ONT"

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
    COLLATE_FASTQ(meta_cram_ch )
    | set { reads_ch }

    if (params.cleanup_intermediate_files_irods_extractor) {
        COLLATE_FASTQ.out.remove_channel.flatten()
                .filter(Path)
                .map { it.delete() }
    }

    emit: reads_ch // tuple val(meta), path(forward_fastq), path(reverse_fastq
}