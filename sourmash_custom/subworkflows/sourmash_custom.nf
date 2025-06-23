#!/usr/bin/env nextflow

/*
========================================================================================
    IMPORT MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULES
//
include { SOURMASH_SKETCH_DB } from '../modules/sourmash.nf'
include { DREP_GENERATE_STB  } from '../modules/drep.nf'

workflow SOURMASH_DATABASE {
    take:
    assembly_manifest_ch        

    main: 

    DREP_GENERATE_STB(params.db_name, assembly_manifest_ch)
    SOURMASH_SKETCH_DB(params.db_name, assembly_manifest_ch)

}

workflow {
    assembly_manifest_ch = Channel.fromPath(params.db_manifest)
    SOURMASH_DATABASE(assembly_manifest_ch)
}