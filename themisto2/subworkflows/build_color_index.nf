include { COLOR_MAPPING                                    } from '../modules/color_mapping.nf'
include { GGCAT as GGCAT_SPECIES                            } from '../modules/ggcat.nf'
include { GGCAT as GGCAT_GROUP                              } from '../modules/ggcat.nf'
include { GGCAT as GGCAT_CANDIDATE                          } from '../modules/ggcat.nf'
include { SBWT_BUILD as SBWT_BUILD_SPECIES; SBWT_CHECK as SBWT_CHECK_SPECIES } from '../modules/sbwt.nf'
include { SBWT_BUILD as SBWT_BUILD_GROUP; SBWT_CHECK as SBWT_CHECK_GROUP } from '../modules/sbwt.nf'
include { SBWT_BUILD as SBWT_BUILD_CANDIDATE; SBWT_CHECK as SBWT_CHECK_CANDIDATE } from '../modules/sbwt.nf'
include { THEMISTO2_BUILD as THEMISTO2_BUILD_SPECIES; THEMISTO2_STATS as THEMISTO2_STATS_SPECIES; THEMISTO2_EXPORT as THEMISTO2_EXPORT_SPECIES } from '../modules/themisto2.nf'
include { THEMISTO2_BUILD as THEMISTO2_BUILD_CANDIDATE; THEMISTO2_STATS as THEMISTO2_STATS_CANDIDATE } from '../modules/themisto2.nf'
include { LINEAGE_SPECIFICITY_FILTER; CANDIDATE_COLOR_LIST  } from '../modules/lineage_specificity_filtering.nf'
include { validate_parameters                               } from '../modules/validate_parameters.nf'

workflow BUILD_COLOR_INDEX {
    take:
    // One pre-paired item per species: [ meta, metadata, assembly_input ]
    //   meta.ID            -- species name; the output folder and the join key
    //                         between species_index and lineage_index
    //   meta.target_groups -- comma-separated lineage labels for this species
    //                         (empty/absent = species-wide only)
    // The caller pairs metadata+assembly and sets per-species target_groups (see
    // lsmd subworkflows/manifest_parse.nf) -- replaces the old combine() cross-
    // multiply and the single global params.target_groups.
    samples_ch

    main:
    validate_parameters()

    color_mapping_input = samples_ch // [meta, metadata, assembly_input]

    // Step 02 - metadata + assemblies -> Themisto colour-file format
    COLOR_MAPPING(color_mapping_input)

    // ============ Index A (species-wide) -- always built ============

    // Step 03 - unitigs from the colour file
    GGCAT_SPECIES(COLOR_MAPPING.out.file_colors)

    // Step 04 - build SBWT, then verify (SBWT_CHECK split out: much lighter than the
    // build; generic, so also reused downstream). Check the .sbwt alone, re-pair .lcs after.
    SBWT_BUILD_SPECIES(GGCAT_SPECIES.out.unitigs)

    SBWT_BUILD_SPECIES.out.index
    | map { meta, sbwt, lcs -> [meta, sbwt] }
    | set { sbwt_only }

    SBWT_CHECK_SPECIES(sbwt_only)

    SBWT_CHECK_SPECIES.out.index
    | join(SBWT_BUILD_SPECIES.out.index.map { meta, sbwt, lcs -> [meta, lcs] })
    | set { checked_index } // tuple(meta, sbwt, lcs)

    // Step 05 - Themisto2 index (colour file + verified SBWT index)
    COLOR_MAPPING.out.file_colors
    | join(checked_index)
    | set { themisto_build_input }

    THEMISTO2_BUILD_SPECIES(themisto_build_input)
    THEMISTO2_STATS_SPECIES(THEMISTO2_BUILD_SPECIES.out.index)

    // Step 06 - export
    THEMISTO2_EXPORT_SPECIES(THEMISTO2_STATS_SPECIES.out.index)

    // ============ Index B (per-lineage) -- only when --target_groups is set ============
    // Same steps 03-06, run as a separate aliased pipeline rather than merging into A's channels.

    // Fan COLOR_MAPPING's per-lineage outputs into one item per lineage, tagged
    // meta = [ID: <group>, species: <species run ID>]. The .species key is how
    // SET_DIFF_CALCULATIONS joins lineage_index against species_index. The glob emits
    // a List when >1 group matched, a bare Path when exactly 1 -- handle both.
    def fan_out_target_group = { meta, files ->
        def file_list = files instanceof List ? files : [files]
        file_list.collect { f -> [[ID: f.getParent().getName(), species: meta.ID], f] }
    }

    COLOR_MAPPING.out.target_group_file_colors
    | flatMap(fan_out_target_group)
    | set { lineage_file_colors }

    GGCAT_GROUP(lineage_file_colors)

    SBWT_BUILD_GROUP(GGCAT_GROUP.out.unitigs)

    SBWT_BUILD_GROUP.out.index
    | map { meta, sbwt, lcs -> [meta, sbwt] }
    | set { lineage_sbwt_only }

    SBWT_CHECK_GROUP(lineage_sbwt_only)

    SBWT_CHECK_GROUP.out.index
    | join(SBWT_BUILD_GROUP.out.index.map { meta, sbwt, lcs -> [meta, lcs] })
    | set { lineage_checked_index } // tuple(meta, sbwt, lcs), meta.species set -- Index B

    // ============ Step 07 - lineage-specificity candidate filtering ============
    // lineage_specificity_filter.py runs ONCE per species over the SPECIES-wide export
    // (not Index B's lineage-scoped export): keeps unitigs that are lineage-core
    // (within_frac) and, when --specificity_max_outside is set, lineage-specific (max
    // presence across any single OTHER lineage). Emits one candidate FASTA per lineage
    // in meta.target_groups.
    //
    // Gated on meta.target_groups: the species export is always non-empty, so a species
    // with no requested lineages must be filtered out here rather than relying on an
    // empty upstream channel (as the Index B branch above does via its optional glob).
    THEMISTO2_EXPORT_SPECIES.out.unitigs
    | join(THEMISTO2_EXPORT_SPECIES.out.color_sets)
    | join(THEMISTO2_EXPORT_SPECIES.out.metadata)
    | join(COLOR_MAPPING.out.label_mapping)
    | filter { meta, unitigs, color_sets, export_metadata, label_mapping -> meta.target_groups }
    | set { specificity_filter_input } // tuple(meta, unitigs, color_sets, export_metadata, label_mapping) -- species-wide

    LINEAGE_SPECIFICITY_FILTER(specificity_filter_input)

    // Fan the per-lineage outputs out the same way fan_out_target_group does above
    // (meta = [ID: lineage, species: species run ID], needed for the SET_DIFF join) --
    // here the lineage ID comes from the filename, not a colour-file subdirectory.
    LINEAGE_SPECIFICITY_FILTER.out.unitigs
    | flatMap { meta, files ->
        def file_list = files instanceof List ? files : [files]
        file_list.collect { f -> [[ID: (f.name - '_candidate_unitigs.fasta'), species: meta.ID], f] }
    }
    // Empty FASTA = nothing cleared the thresholds -- skip the rebuild for that lineage.
    // log.warn so it shows in the main pipeline log, not just the task work-dir.
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
    | set { candidate_checked_index } // tuple(meta, sbwt, lcs), meta.species set -- E

    // QC gate only -- confirms the rebuild is structurally sound. No export: the
    // candidate index only feeds sbwt difference, which never reads exported unitigs.
    CANDIDATE_COLOR_LIST.out.file_colors
    | join(candidate_checked_index)
    | set { candidate_themisto_build_input }

    THEMISTO2_BUILD_CANDIDATE(candidate_themisto_build_input)
    THEMISTO2_STATS_CANDIDATE(THEMISTO2_BUILD_CANDIDATE.out.index)

    emit:
    // Public contract = only the SBWT indexes set_diff_calculations.nf needs. The species
    // Themisto2 export still runs above -- LINEAGE_SPECIFICITY_FILTER (step 07) reads it.
    sbwt_index       = checked_index           // tuple(meta, sbwt, lcs) -- species-wide (A)
    lineage_index    = lineage_checked_index   // tuple(meta, sbwt, lcs), meta.species set -- B
    candidate_index  = candidate_checked_index // tuple(meta, sbwt, lcs), meta.species set -- E
}
