include { COLLATE_CRAM; FASTQ_FROM_COLLATED_BAM } from '../modules/irods/samtools.nf'
include { BATON } from '../modules/irods/baton.nf'
include { JSON_PREP; JSON_PARSE } from '../modules/irods/jq.nf'
include { RETRIEVE_CRAM } from '../modules/irods/retrieve.nf'

workflow IRODS_QUERY {
        take:
        input_irods_ch //tuple study, runid

        main:
        JSON_PREP(input_irods_ch)
        | BATON
        | JSON_PARSE

        JSON_PARSE.out.paths.splitText().map{ cram_path ->
        def meta = [:]
        meta.ID = cram_path.split("/")[-1].split(".cram")[0]
        [ meta, cram_path ]
        }.set{ meta_cram_ch }

        emit:
        meta_cram_ch

}

workflow IRODS_EXTRACT {
    
    take:
    input_irods_ch //tuple study, runid

    main:
    IRODS_QUERY(input_irods_ch).set{ meta_cram_ch }

    Channel.fromPath("${params.outdir}/*#*/raw_fastq/*_1.fastq.gz").map{ raw_fastq_path ->
        ID = raw_fastq_path.simpleName.split("_1")[0]
    }.ifEmpty("fresh_run").set{ existing_id }

    meta_cram_ch.combine( existing_id | collect | map{ [it] })
    | filter { meta, cram_path, existing -> !(meta.ID in existing)}
    | view
    | map { it[0,1] }
    | set{ do_not_exist }

    RETRIEVE_CRAM(do_not_exist)
    | COLLATE_CRAM
    | FASTQ_FROM_COLLATED_BAM

    FASTQ_FROM_COLLATED_BAM.out.remove_channel.flatten()
            .filter(Path)
            .map { it.delete() }

    emit:
    reads_ch = FASTQ_FROM_COLLATED_BAM.out.fastq_channel
}