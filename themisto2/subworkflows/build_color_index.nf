// Species-wide colour index build ONLY -- GGCAT unitigs -> SBWT -> Themisto2
// build/stats/export. No filtering of any kind happens here any more: lineage-specificity
// filtering and the candidate-index rebuild moved to marker_filtering.nf, which also owns
// the ATB cross-species check that used to be the sbwt bg_excl/markers set-diff. This split
// keeps a species-index rebuild's cache untouched by edits to filtering params/thresholds,
// and vice versa.
include { COLOR_MAPPING                                    } from '../modules/color_mapping.nf'
include { GGCAT as GGCAT_SPECIES                            } from '../modules/ggcat.nf'
include { SBWT_BUILD as SBWT_BUILD_SPECIES; SBWT_CHECK as SBWT_CHECK_SPECIES } from '../modules/sbwt.nf'
include { THEMISTO2_BUILD as THEMISTO2_BUILD_SPECIES; THEMISTO2_STATS as THEMISTO2_STATS_SPECIES; THEMISTO2_EXPORT as THEMISTO2_EXPORT_SPECIES } from '../modules/themisto2.nf'
include { validate_parameters                               } from '../modules/validate_parameters.nf'
include { CHECKPOINT_COUNT                                  } from '../modules/checkpoint_count.nf'

workflow BUILD_COLOR_INDEX {
    take:
    samples_ch          // [ [ID: species], metadata, assembly_input ]

    main:
    validate_parameters()

    color_mapping_input = samples_ch // [meta, metadata, assembly_input]

    // Metadata + assemblies -> Themisto colour-file format
    COLOR_MAPPING(color_mapping_input)

    // ============ Species-wide index -- always built ============

    // Unitigs from the colour file
    GGCAT_SPECIES(COLOR_MAPPING.out.file_colors)

    // Build SBWT, then verify (SBWT_CHECK split out: much lighter than the
    // build; generic, so also reused downstream). Check the .sbwt alone, re-pair .lcs after.
    SBWT_BUILD_SPECIES(GGCAT_SPECIES.out.unitigs)

    SBWT_BUILD_SPECIES.out.index
    | map { meta, sbwt, lcs -> [meta, sbwt] }
    | set { sbwt_only }

    SBWT_CHECK_SPECIES(sbwt_only)

    SBWT_CHECK_SPECIES.out.index
    | join(SBWT_BUILD_SPECIES.out.index.map { meta, sbwt, lcs -> [meta, lcs] })
    | set { checked_index } // tuple(meta, sbwt, lcs)

    // Themisto2 index (colour file + verified SBWT index)
    COLOR_MAPPING.out.file_colors
    | join(checked_index)
    | set { themisto_build_input }

    THEMISTO2_BUILD_SPECIES(themisto_build_input)
    THEMISTO2_STATS_SPECIES(THEMISTO2_BUILD_SPECIES.out.index)

    // Export
    THEMISTO2_EXPORT_SPECIES(THEMISTO2_STATS_SPECIES.out.index)

    THEMISTO2_EXPORT_SPECIES.out.unitigs
    | join(THEMISTO2_EXPORT_SPECIES.out.color_sets)
    | join(THEMISTO2_EXPORT_SPECIES.out.metadata)
    | join(COLOR_MAPPING.out.label_mapping)
    | set { species_export } // tuple(meta, unitigs, color_sets, export_metadata, label_mapping) -- species-wide

    // ============ Checkpoint counts ============
    // Side-channel only: tap each key stage's output, never joined back in. `order` keys
    // leave room (gaps of 10) for marker_filtering.nf's stages to continue the funnel.
    Channel.empty()
    | mix( COLOR_MAPPING.out.file_colors.map          { meta, f -> [meta, 10, 'species_colorfile', 'colorfile', f] } )
    | mix( GGCAT_SPECIES.out.unitigs.map              { meta, f -> [meta, 20, 'species_ggcat_unitigs', 'fasta', f] } )
    | mix( THEMISTO2_BUILD_SPECIES.out.index.map      { meta, f -> [meta, 30, 'species_themisto_index', 'themisto', f] } )
    | mix( THEMISTO2_EXPORT_SPECIES.out.unitigs.map   { meta, f -> [meta, 40, 'species_export_unitigs', 'fasta', f] } )
    | set { checkpoint_inputs }

    CHECKPOINT_COUNT(checkpoint_inputs)

    emit:
    sbwt_index      = checked_index  // tuple(meta, sbwt, lcs) -- species-wide
    species_export  = species_export // tuple(meta, unitigs, color_sets, export_metadata, label_mapping) -- marker_filtering.nf's LINEAGE_SPECIFICITY_FILTER input
    checkpoints     = CHECKPOINT_COUNT.out.row // tuple(meta, row_tsv) -- per-stage count rows
}
