include { IRODS_MANIFEST_PARSE     } from './subworkflows/irods_manifest_parse.nf'
include { INPUT_CHECK              } from './subworkflows/input_check.nf'
include { IRODS_CLI; COMBINE_IRODS } from './subworkflows/combined_input.nf'
include { ENA_DOWNLOAD             } from './subworkflows/ena_input.nf'
include { IRODS_EXTRACTOR          } from '../irods_extractor/subworkflows/irods.nf'

include { validate_parameters      } from './modules/validate_parameters'

workflow MIXED_INPUT {
    /*
    This workflow is a handler to call only what is needed from the IRODS workflows
    replaces the common stucture of 
    COMBINE_IRODS
    | IRODS_EXTRACTOR
    | COMBINE_READS

    from the start of IRODS_PIPELINES and allows for smarter and validated use of the params given
    */

    // No inputs
    main:
    def active_workflows = validate_parameters() //ensure all potential inputs are parsed and validated

    if ('ENA' in active_workflows) {
        Channel.fromPath(params.manifest_ena)
        | ENA_DOWNLOAD
        | set { reads_from_ena_ch }
    } else {
        Channel.of("none")
        | set { reads_from_ena_ch }
    }

    if ('IRODS' in active_workflows) {
        COMBINE_IRODS
        | IRODS_EXTRACTOR
        | set { reads_from_irods_ch }
    } else {
        Channel.of("none")
        | set { reads_from_irods_ch }
    }

    if ('READS_MANIFEST' in active_workflows ) {
        def manifestToUse = params.manifest_of_reads ? params.manifest_of_reads : params.manifest

        input_reads_ch = file(manifestToUse)

        INPUT_CHECK(input_reads_ch)
        | set { reads_from_local_ch }
    } else {
        Channel.of("none")
        | set { reads_from_local_ch }
    }

    reads_from_irods_ch
    | mix(reads_from_local_ch)
    | mix(reads_from_ena_ch)
    | filter { it != "none"}
    | set { all_reads_ready_ch }

    emit:
    all_reads_ready_ch //channel of [meta, R1, R2] taken from a mixture of IRODS + Given manifest or either
}