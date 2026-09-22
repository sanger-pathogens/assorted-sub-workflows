include { COLOR_MAPPING                                    } from '../modules/color_mapping.nf'
include { GGCAT as GGCAT_SPECIES                            } from '../modules/ggcat.nf'
include { SBWT_BUILD as SBWT_BUILD_SPECIES; SBWT_CHECK as SBWT_CHECK_SPECIES } from '../modules/sbwt.nf'
include { THEMISTO2_BUILD as THEMISTO2_BUILD_SPECIES; THEMISTO2_STATS as THEMISTO2_STATS_SPECIES; THEMISTO2_EXPORT as THEMISTO2_EXPORT_SPECIES } from '../modules/themisto2.nf'
include { validate_parameters                               } from '../modules/validate_parameters.nf'
include { CHECKPOINT_FASTA; CHECKPOINT_THEMISTO             } from '../modules/checkpoint.nf'

workflow BUILD_COLOR_INDEX {
    take:
    samples_ch

    main:
    validate_parameters()

    color_mapping_input = samples_ch

    COLOR_MAPPING(color_mapping_input)

    GGCAT_SPECIES(COLOR_MAPPING.out.file_colors)

    SBWT_BUILD_SPECIES(GGCAT_SPECIES.out.unitigs)

    SBWT_BUILD_SPECIES.out.index
    | map { meta, sbwt, lcs -> [meta, sbwt] }
    | set { sbwt_only }

    SBWT_CHECK_SPECIES(sbwt_only)

    SBWT_CHECK_SPECIES.out.index
    | join(SBWT_BUILD_SPECIES.out.index.map { meta, sbwt, lcs -> [meta, lcs] })
    | set { checked_index }

    COLOR_MAPPING.out.file_colors
    | join(checked_index)
    | set { themisto_build_input }

    THEMISTO2_BUILD_SPECIES(themisto_build_input)
    THEMISTO2_STATS_SPECIES(THEMISTO2_BUILD_SPECIES.out.index)

    THEMISTO2_EXPORT_SPECIES(THEMISTO2_STATS_SPECIES.out.index)

    THEMISTO2_EXPORT_SPECIES.out.unitigs
    | join(THEMISTO2_EXPORT_SPECIES.out.color_sets)
    | join(THEMISTO2_EXPORT_SPECIES.out.metadata)
    | join(COLOR_MAPPING.out.label_mapping)
    | set { species_export }

    Channel.empty()
    | mix( COLOR_MAPPING.out.file_colors.map          { meta, f -> [meta, 10, 'species_colorfile', 'colorfile', f] } )
    | mix( GGCAT_SPECIES.out.unitigs.map              { meta, f -> [meta, 20, 'species_ggcat_unitigs', 'fasta', f] } )
    | mix( THEMISTO2_EXPORT_SPECIES.out.unitigs.map   { meta, f -> [meta, 40, 'species_export_unitigs', 'fasta', f] } )
    | set { fasta_checkpoint_inputs }

    THEMISTO2_BUILD_SPECIES.out.index.map { meta, f -> [meta, 30, 'species_themisto_index', 'themisto', f] }
    | set { themisto_checkpoint_inputs }

    CHECKPOINT_FASTA(fasta_checkpoint_inputs)
    CHECKPOINT_THEMISTO(themisto_checkpoint_inputs)

    emit:
    sbwt_index      = checked_index
    species_export  = species_export
    checkpoints     = CHECKPOINT_FASTA.out.row.mix(CHECKPOINT_THEMISTO.out.row)
}
