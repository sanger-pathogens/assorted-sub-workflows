include { COLOR_MAPPING                                    } from '../modules/color_mapping.nf'
include { GGCAT as GGCAT_SPECIES                            } from '../modules/ggcat.nf'
include { GGCAT as GGCAT_CANDIDATE                          } from '../modules/ggcat.nf'
include { SBWT_BUILD as SBWT_BUILD_SPECIES; SBWT_CHECK as SBWT_CHECK_SPECIES } from '../modules/sbwt.nf'
include { SBWT_BUILD as SBWT_BUILD_CANDIDATE; SBWT_CHECK as SBWT_CHECK_CANDIDATE } from '../modules/sbwt.nf'
include { THEMISTO2_BUILD as THEMISTO2_BUILD_SPECIES; THEMISTO2_STATS as THEMISTO2_STATS_SPECIES; THEMISTO2_EXPORT as THEMISTO2_EXPORT_SPECIES } from '../modules/themisto2.nf'
include { THEMISTO2_BUILD as THEMISTO2_BUILD_CANDIDATE; THEMISTO2_STATS as THEMISTO2_STATS_CANDIDATE } from '../modules/themisto2.nf'
include { LINEAGE_SPECIFICITY_FILTER; CANDIDATE_COLOR_LIST  } from '../modules/lineage_specificity_filtering.nf'
include { validate_parameters                               } from '../modules/validate_parameters.nf'
include { CHECKPOINT_COUNT                                  } from '../modules/checkpoint_count.nf'

workflow BUILD_COLOR_INDEX {
    take:
    samples_ch

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

    // ============ Lineage-specificity candidate filtering ============
    // lineage_specificity_filter.py runs ONCE per species over the SPECIES-wide export:
    // keeps unitigs that are lineage-core (within_frac) and, when
    // --specificity_max_outside is set, lineage-specific (max presence across any single
    // OTHER lineage). Emits one candidate FASTA per targeted lineage.
    //
    // Which lineages are targeted comes from meta.target_groups: a non-empty value
    // targets exactly those lineages; blank/absent targets every lineage with
    // >= candidate_min_genome_count genomes (excluding 'unclassified'). Runs for every
    // species -- a species with no usable lineage just emits nothing (outputs optional).
    THEMISTO2_EXPORT_SPECIES.out.unitigs
    | join(THEMISTO2_EXPORT_SPECIES.out.color_sets)
    | join(THEMISTO2_EXPORT_SPECIES.out.metadata)
    | join(COLOR_MAPPING.out.label_mapping)
    | set { specificity_filter_input } // tuple(meta, unitigs, color_sets, export_metadata, label_mapping) -- species-wide

    LINEAGE_SPECIFICITY_FILTER(specificity_filter_input)

    // One task emits one candidate FASTA per requested lineage. Fan them into one
    // item each, meta = [ID: lineage, species: species run ID] (meta.species is the
    // SET_DIFF join key) -- the lineage ID comes from the filename.
    LINEAGE_SPECIFICITY_FILTER.out.unitigs
    | flatMap { meta, files ->
        def file_list = files instanceof List ? files : [files]
        file_list.collect { f -> [[ID: (f.name - '_candidate_unitigs.fasta'), species: meta.ID], f] }
    }
    | set { candidate_fasta_per_lineage } // tuple(meta, fasta) -- one per lineage, may be empty

    // Empty FASTA = nothing cleared the thresholds -- skip the rebuild for that lineage.
    // log.warn so it shows in the main pipeline log, not just the task work-dir.
    candidate_fasta_per_lineage
    | filter { meta, fasta ->
        if (fasta.size() == 0) {
            log.warn("No candidate unitigs survived filtering for lineage '${meta.ID}' -- skipping candidate_index rebuild.")
            return false
        }
        true
    }
    // stage: 'candidate' -- meta.ID/species stay the lineage's own (needed for the
    // SET_DIFF join); the stage key is what disambiguates this rebuild's publishDir.
    | map { meta, fasta -> [meta + [stage: 'candidate'], fasta] }
    | set { candidate_fasta_nonempty }

    // Rebuild the candidate FASTA into a real index. CANDIDATE_COLOR_LIST wraps it as a
    // one-line colour-list file first -- GGCAT/THEMISTO2_BUILD expect that, not a raw FASTA.
    CANDIDATE_COLOR_LIST(candidate_fasta_nonempty)

    GGCAT_CANDIDATE(CANDIDATE_COLOR_LIST.out.file_colors)
    SBWT_BUILD_CANDIDATE(GGCAT_CANDIDATE.out.unitigs)

    SBWT_BUILD_CANDIDATE.out.index
    | map { meta, sbwt, lcs -> [meta, sbwt] }
    | set { candidate_sbwt_only }

    SBWT_CHECK_CANDIDATE(candidate_sbwt_only)

    SBWT_CHECK_CANDIDATE.out.index
    | join(SBWT_BUILD_CANDIDATE.out.index.map { meta, sbwt, lcs -> [meta, lcs] })
    | set { candidate_checked_index } // tuple(meta, sbwt, lcs), meta.species set -- candidate index

    // QC gate only -- confirms the rebuild is structurally sound. No export: the
    // candidate index only feeds sbwt difference, which never reads exported unitigs.
    CANDIDATE_COLOR_LIST.out.file_colors
    | join(candidate_checked_index)
    | set { candidate_themisto_build_input }

    THEMISTO2_BUILD_CANDIDATE(candidate_themisto_build_input)
    THEMISTO2_STATS_CANDIDATE(THEMISTO2_BUILD_CANDIDATE.out.index)

    // ============ Checkpoint counts ============
    // Side-channel only: tap each key stage's output, never joined back in. Each row
    // carries an ascending `order` key so pipeline_counts.tsv reads top-down as the
    // funnel (species-wide build -> candidate filter -> candidate build). Gaps of 10
    // leave room for the caller's own later stages.
    Channel.empty()
    | mix( COLOR_MAPPING.out.file_colors.map          { meta, f -> [meta, 10, 'species_colorfile', 'colorfile', f] } )
    | mix( GGCAT_SPECIES.out.unitigs.map              { meta, f -> [meta, 20, 'species_ggcat_unitigs', 'fasta', f] } )
    | mix( THEMISTO2_BUILD_SPECIES.out.index.map      { meta, f -> [meta, 30, 'species_themisto_index', 'themisto', f] } )
    | mix( THEMISTO2_EXPORT_SPECIES.out.unitigs.map   { meta, f -> [meta, 40, 'species_export_unitigs', 'fasta', f] } )
    | mix( candidate_fasta_per_lineage.map            { meta, f -> [meta, 50, 'candidate_specificity_filter', 'fasta', f] } )
    | mix( THEMISTO2_BUILD_CANDIDATE.out.index.map    { meta, f -> [meta, 60, 'candidate_rebuilt_index', 'themisto', f] } )
    | set { checkpoint_inputs }

    CHECKPOINT_COUNT(checkpoint_inputs)

    emit:
    // Public contract = only the SBWT indexes set_diff_calculations.nf needs. The species
    // Themisto2 export still runs above -- LINEAGE_SPECIFICITY_FILTER reads it.
    sbwt_index       = checked_index           // tuple(meta, sbwt, lcs) -- species-wide
    candidate_index  = candidate_checked_index // tuple(meta, sbwt, lcs), meta.species set -- candidate index
    checkpoints      = CHECKPOINT_COUNT.out.row // tuple(meta, row_tsv) -- per-stage count rows
}
