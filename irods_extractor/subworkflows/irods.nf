include { COLLATE_FASTQ; COMBINE_FASTQ  } from '../modules/samtools.nf'
include { BATON                         } from '../modules/baton.nf'
include { JSON_PREP                     } from '../modules/jq.nf'
include { METADATA as METADATA_QUERIED  } from '../modules/metadata_save.nf'
include { CRAM_EXTRACT                  } from './extraction_methods/illumina_extract.nf'

include { ILLUMINA_PARSE               } from './extraction_methods/illumina_extract.nf'
include { ONT_PARSE                    } from './extraction_methods/ONT_extract.nf'

def IRODS_ERROR_MSG = """
    Error: IRODS search returned no data!
    - Please ensure you are logged on with `iinit` and re-run the
    pipeline without `-resume`.
    - If you are logged in, check the process IRODS_QUERY:BATON work
    directory for permission errors and ensure the study exists.
    - ensure your --read_type (${params.read_type}) is correct
"""

workflow IRODS_QUERY {
        take:
        input_irods_ch // map

        main:
        if (params.fastsearch) {
            input_irods_ch.collect()
            | set { search_ch }
        } else {
            input_irods_ch
            | set { search_ch }
        }

        JSON_PREP(search_ch)
        | BATON
        | splitJson(path: "result.multiple")
        | branch { meta ->
            def object = meta?.data_object
            extract: (object && (object.endsWith('.cram') || object.endsWith('.bam')))
            collection: (!object)
        }
        | set { data_type }
        
        ILLUMINA_PARSE(data_type.extract)
        | set{ illumina_ch }

        ONT_PARSE(data_type.collection)
        | set{ ont_ch }

        illumina_ch
        | mix(ont_ch)
        | filter { it && it instanceof Map && !it.isEmpty() } //make sure it is a non-empty map at this stage
        | ifEmpty { error(IRODS_ERROR_MSG) }
        | set { meta_ch }

        if (params.save_metadata) {
            meta_ch
            | collectFile() { map -> [ "lane_metadata.txt", map.collect{it}.join(', ') + '\n' ] }
            | set{ metadata_only }

            METADATA_QUERIED(metadata_only, "irods_queried")
        }

        emit:
        illumina_ch
        ont_ch
}

workflow IRODS_EXTRACTOR {
    // THIS IS THE WITHIN PIPELINE IRODS EXTRACTOR
    take:
    input_irods_ch // map

    main:
    //todo remove or make new method for ONT
    if (!params.illumina) {
        log.error("Only Illumina reads are supported in this pipeline")
    }
    //expects short reads currently.
    IRODS_QUERY(input_irods_ch)
    | CRAM_EXTRACT

    emit:
    reads_ch = CRAM_EXTRACT.out.reads_ch // tuple val(meta), path(forward_fastq), path(reverse_fastq)
}