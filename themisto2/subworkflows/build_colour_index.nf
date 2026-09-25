include { COLOUR_MAPPING                                   } from '../modules/colour_mapping.nf'
include { GGCAT as GGCAT_SPECIES                            } from '../modules/ggcat.nf'
include { SBWT_BUILD as SBWT_BUILD_SPECIES; SBWT_CHECK as SBWT_CHECK_SPECIES; SBWT_DUMP_UNITIGS as SBWT_DUMP_UNITIGS_SPECIES } from '../modules/sbwt.nf'
include { THEMISTO2_BUILD as THEMISTO2_BUILD_SPECIES; THEMISTO2_STATS as THEMISTO2_STATS_SPECIES; THEMISTO2_EXPORT as THEMISTO2_EXPORT_SPECIES } from '../modules/themisto2.nf'
include { validate_parameters                               } from '../modules/validate_parameters.nf'
include { CHECKPOINT_FASTA; CHECKPOINT_THEMISTO             } from '../modules/checkpoint.nf'

workflow BUILD_COLOUR_INDEX {
    take:
    samples_ch

    main:
    validate_parameters()

    colour_mapping_input = samples_ch

    COLOUR_MAPPING(colour_mapping_input)

    GGCAT_SPECIES(COLOUR_MAPPING.out.file_colours)

    SBWT_BUILD_SPECIES(GGCAT_SPECIES.out.unitigs)

    SBWT_BUILD_SPECIES.out.index
    | map { meta, sbwt, lcs -> [meta, sbwt] }
    | set { sbwt_only }

    SBWT_CHECK_SPECIES(sbwt_only)

    // Pass or fail goes to .nextflow.log; a failed index goes no further.
    SBWT_CHECK_SPECIES.out.index
    | filter { meta, sbwt, result ->
        if (result == 'PASS') {
            log.info("SBWT check passed: species index ${meta.ID}")
            return true
        }
        log.warn("SBWT check FAILED: species index ${meta.ID} -- skipping everything downstream for this species.")
        false
    }
    | map { meta, sbwt, result -> [meta, sbwt] }
    | set { species_sbwt_checked }

    // Checkpoint only: both strands of every unitig (the index is built with -r).
    SBWT_DUMP_UNITIGS_SPECIES(species_sbwt_checked)

    species_sbwt_checked
    | join(SBWT_BUILD_SPECIES.out.index.map { meta, sbwt, lcs -> [meta, lcs] })
    | set { checked_index }

    COLOUR_MAPPING.out.file_colours
    | join(checked_index)
    | set { themisto_build_input }

    THEMISTO2_BUILD_SPECIES(themisto_build_input)
    THEMISTO2_STATS_SPECIES(THEMISTO2_BUILD_SPECIES.out.index)

    THEMISTO2_EXPORT_SPECIES(THEMISTO2_STATS_SPECIES.out.index)

    THEMISTO2_EXPORT_SPECIES.out.unitigs
    | join(THEMISTO2_EXPORT_SPECIES.out.colour_sets)
    | join(THEMISTO2_EXPORT_SPECIES.out.metadata)
    | join(COLOUR_MAPPING.out.label_mapping)
    | set { species_export }

    // Stage names and their order in the report live in checkpoint.nf (checkpoint_steps()).
    Channel.empty()
    | mix( COLOUR_MAPPING.out.file_colours.map          { meta, f -> [meta, 'species_colourfile', 'colourfile', f] } )
    | mix( GGCAT_SPECIES.out.unitigs.map                { meta, f -> [meta, 'species_ggcat_unitigs', 'fasta', f] } )
    | mix( SBWT_DUMP_UNITIGS_SPECIES.out.unitigs.map    { meta, f -> [meta, 'species_sbwt_unitigs', 'fasta', f] } )
    | mix( THEMISTO2_EXPORT_SPECIES.out.unitigs.map     { meta, f -> [meta, 'species_export_unitigs', 'fasta', f] } )
    | set { fasta_checkpoint_inputs }

    THEMISTO2_STATS_SPECIES.out.stats.map { meta, f -> [meta, 'species_themisto_index', f] }
    | set { themisto_checkpoint_inputs }

    CHECKPOINT_FASTA(fasta_checkpoint_inputs)
    CHECKPOINT_THEMISTO(themisto_checkpoint_inputs)

    emit:
    sbwt_index      = checked_index
    species_export  = species_export
    checkpoints     = CHECKPOINT_FASTA.out.row.mix(CHECKPOINT_THEMISTO.out.row)
}
