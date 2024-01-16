include { IRODS_MANIFEST_PARSE } from '../../irods_extractor/subworkflows/irods_manifest_parse.nf'

//
// SUBWORKFLOW: Read in study, run, etc. parameters and pull data from iRODS
//
workflow IRODS_CLI{
    main:
    param_input = Channel.of(["${params.studyid}", "${params.runid}", "${params.laneid}", "${params.plexid}"])
    
    param_input.map{ studyid, runid, laneid, plexid ->
        meta = [:]
        if (studyid > 0) {meta.studyid = studyid}
        if (runid > 0) {meta.runid = runid}
        if (laneid > 0 ) {meta.laneid = laneid}
        if (plexid > 0 ) {meta.plexid = plexid}
        meta
    }.set{ input_irods_from_opt_ch } 

    emit:
    input_irods_from_opt_ch
}

workflow COMBINE_IRODS{
    main:
    // take iRODS dataset specification from CLI options
    if (params.studyid > 0) {
        IRODS_CLI()
    } else {
        Channel.of("none").set{ input_irods_from_opt_ch }
    }

    // take iRODS dataset specification from manifest of lanes
    if (params.manifest_of_lanes) {
        IRODS_MANIFEST_PARSE(params.manifest_of_lanes)
        | set{ input_irods_from_man_ch }
    } else {
        Channel.of("none").set{ input_irods_from_man_ch}
    }

    // combine iRODS specs input channels
    input_irods_from_opt_ch.mix(input_irods_from_man_ch).set{ input_irods_ch }

    emit:
    input_irods_ch
}

workflow COMBINE_READS{
    take:
    reads_ch // [meta, read_1, read_2] as from IRODS_EXTRACTOR

    // change channel structure to match that from INPUT_CHECK
    reads_ch.map{ meta, read_1, read_2 ->
        meta_new = [:]
        meta_new.id = meta.ID
        reads = [read_1, read_2]
        [ meta_new, reads ]
    }.set{ irods_ready_to_map_ch }

    // combine reads input channels
    irods_ready_to_map_ch.mix(ch_reads_from_manifest.filter{ it != "none"}).set{ all_reads_ready_to_map_ch }

    emit:
    all_reads_ready_to_map_ch

}