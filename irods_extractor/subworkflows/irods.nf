include { COLLATE_FASTQ; COMBINE_FASTQ  } from '../modules/samtools.nf'
include { BATON                         } from '../modules/baton.nf'
include { JSON_PREP; JSON_PARSE         } from '../modules/jq.nf'
include { METADATA as METADATA_QUERIED  } from '../modules/metadata_save.nf'
include { CRAM_EXTRACT                  } from './extraction_methods/illumina_extract.nf'

include { ILLUMINA_PARSE                } from './extraction_methods/illumina_extract.nf'
include { ONT_PARSE                     } from './extraction_methods/ONT_extract.nf'
include { FILTER_EXISTING_OUTPUTS       } from './skip_downloaded.nf'

workflow IRODS_QUERY {
        take:
        input_irods_ch // map

        main:
        JSON_PREP(input_irods_ch)
        | BATON
        | JSON_PARSE

        switch (params.read_type.toLowerCase()) {
            case "illumina":
                ILLUMINA_PARSE(JSON_PARSE.out.json_file)
                | set{ meta_file_ch }

                break

            case "ont":
                ONT_PARSE(JSON_PARSE.out.json_file)
                | set{ meta_file_ch }

                break

            default:
                log.error("input --read_type: ${params.read_type} was not one of Illumina|ont")
        }
        

        if (params.save_metadata) {
            meta_file_ch.map{metadata_map, path -> metadata_map}
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
    //todo remove or make new method for ONT
    if (params.read_type.toLowerCase() != "illumina") {
        log.error("Only Illumina reads are supported in this pipeline")
    }
    //expects short reads currently.
    IRODS_QUERY(input_irods_ch)
    | CRAM_EXTRACT

    emit:
    reads_ch = CRAM_EXTRACT.out.reads_ch // tuple val(meta), path(forward_fastq), path(reverse_fastq)
}