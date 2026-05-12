//
// Check local folder of reads and generate samplesheet from it and derive read channels
//
include { MANIFEST_GENERATOR       } from './modules/manifest_generator.nf'
include { FILTER_EXISTING_OUTPUTS       } from './../../../pipeline_chaining/subworkflows/skip_downloaded.nf'

workflow MANIFEST_FROM_DIR {
    take:
    dir_path

    main:
        MANIFEST_GENERATOR(dir_path)

        ch_manifest = MANIFEST_GENERATOR.out.ch_manifest_from_dir

        ch_reads = ch_manifest.splitCsv(header: true)
            .map { row -> tuple(
                    [ ID : row.ID ],
                    file(row.R1),
                    file(row.R2)
                )
            } 
        | set { reads_from_local_dir_ch }

        if (params.only_new_input){
            FILTER_EXISTING_OUTPUTS(reads_from_local_dir_ch)
            FILTER_EXISTING_OUTPUTS.out.do_not_exist
            .set { reads_from_local_dir_to_process }
        } else{
            reads_from_local_dir_to_process = reads_from_local_dir_ch
        }