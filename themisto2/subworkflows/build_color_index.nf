include { COLOR_MAPPING                                    } from '../modules/color_mapping.nf'
include { GGCAT as GGCAT_SPECIES                            } from '../modules/ggcat.nf'
include { GGCAT as GGCAT_GROUP                              } from '../modules/ggcat.nf'
include { GGCAT as GGCAT_CANDIDATE                          } from '../modules/ggcat.nf'
include { SBWT_BUILD as SBWT_BUILD_SPECIES; SBWT_CHECK as SBWT_CHECK_SPECIES } from '../modules/sbwt.nf'
include { SBWT_BUILD as SBWT_BUILD_GROUP; SBWT_CHECK as SBWT_CHECK_GROUP } from '../modules/sbwt.nf'
include { SBWT_BUILD as SBWT_BUILD_CANDIDATE; SBWT_CHECK as SBWT_CHECK_CANDIDATE } from '../modules/sbwt.nf'
include { THEMISTO_BUILD as THEMISTO_BUILD_SPECIES; THEMISTO_STATS as THEMISTO_STATS_SPECIES; THEMISTO_EXPORT as THEMISTO_EXPORT_SPECIES } from '../modules/themisto2.nf'
include { THEMISTO_BUILD as THEMISTO_BUILD_GROUP; THEMISTO_STATS as THEMISTO_STATS_GROUP; THEMISTO_EXPORT as THEMISTO_EXPORT_GROUP } from '../modules/themisto2.nf'
include { THEMISTO_BUILD as THEMISTO_BUILD_CANDIDATE; THEMISTO_STATS as THEMISTO_STATS_CANDIDATE } from '../modules/themisto2.nf'
include { CANDIDATE_FILTER; CANDIDATE_COLOR_LIST            } from '../modules/lineage_index_filtering.nf'
include { validate_parameters                               } from '../modules/validate_parameters.nf'

workflow BUILD_COLOR_INDEX {
    take:
    /// Input files: metadata and assemblies
    metadata_ch
    assembly_ch

    main:
    validate_parameters()


    // TODO: only works for a single run -- combine() cross-multiplies if there's ever
    // >1 metadata/assembly. Switch to a manifest channel of pre-paired tuples for multi-run.
    //
    // TODO (multi-species, separate branch): the natural shape is a samplesheet
    // (species_id, metadata, assembly, target_groups) parsed into a manifest channel
    // -- fixes the combine() issue above AND makes target_groups per-species instead
    // of a single global params.target_groups (V. cholerae vs S. pneumoniae would
    // have entirely different lineage nomenclatures). Also needs: publishDir paths
    // (ggcat.nf/sbwt.nf/themisto2.nf) keyed on meta.ID alone right now -- two species
    // with a coincidentally-identical lineage/group name would collide on disk (even
    // though channel join()s would still be correct, since they key on the whole meta
    // map including .species). See also the temp_dir TODO in ggcat.nf/sbwt.nf, which
    // this makes more likely to actually bite.

    // Tag metadata with an ID and pair with assembly input (plumbing only -- COLOR_MAPPING
    // needs one tuple, not two separate inputs).
    metadata_ch
    | map { metadata -> [["ID": metadata.baseName], metadata] }
    | combine(assembly_ch)
    | set { color_mapping_input }

    // Step 02 - map sample metadata + assemblies to Themisto's colour-file format
    COLOR_MAPPING(color_mapping_input)


    // ===================================================================
    // Index A (species-wide) -- always built. Unchanged from before Index B existed,
    // other than GGCAT/SBWT_BUILD/SBWT_CHECK/THEMISTO_* now being explicitly aliased
    // to *_SPECIES (was bare/unaliased) so every stage below reads the same way.
    // ===================================================================

    // Step 03 - build unitigs from the colour file
    GGCAT_SPECIES(COLOR_MAPPING.out.file_colors)

    // Step 04 - build the SBWT index from the unitigs, then verify it loads correctly
    // (verification split into its own process -- much lighter workload than the build,
    // no reason to hold a cpu_16/mem_32 reservation just to run a quick check).
    // SBWT_CHECK is generic (also reused by set_diff_calculations.nf) and never
    // touches .lcs, so check the .sbwt alone and re-pair it with its (untouched,
    // already-produced-by-build) .lcs afterward for the Themisto2 build below.
    SBWT_BUILD_SPECIES(GGCAT_SPECIES.out.unitigs)

    SBWT_BUILD_SPECIES.out.index
    | map { meta, sbwt, lcs -> [meta, sbwt] }
    | set { sbwt_only }

    SBWT_CHECK_SPECIES(sbwt_only)

    SBWT_CHECK_SPECIES.out.index
    | join(SBWT_BUILD_SPECIES.out.index.map { meta, sbwt, lcs -> [meta, lcs] })
    | set { checked_index } // tuple(meta, sbwt, lcs)

    // Step 05 - build the Themisto2 index (needs the colour file *and* the verified SBWT index)
    COLOR_MAPPING.out.file_colors
    | join(checked_index)
    | set { themisto_build_input }

    THEMISTO_BUILD_SPECIES(themisto_build_input)
    THEMISTO_STATS_SPECIES(THEMISTO_BUILD_SPECIES.out.index)

    // Step 06 - export the index
    THEMISTO_EXPORT_SPECIES(THEMISTO_STATS_SPECIES.out.index)


    // ===================================================================
    // Index B (per-lineage) -- only present when --target_groups is set. Same
    // steps 03-06 as Index A above, run as a separate parallel pipeline (aliased
    // GGCAT_GROUP/SBWT_*_GROUP/THEMISTO_*_GROUP) rather than merging into Index A's
    // channels and splitting the result back apart -- keeps Index A's own code
    // above untouched, and matches the same aliasing pattern already used for the
    // candidate rebuild further down.
    // ===================================================================

    // Fan out COLOR_MAPPING's per-lineage outputs (index_target_group/<group>/...,
    // one process call can match >1 group) into one channel item per lineage, tagged
    // meta = [ID: <group>, species: <species run's meta.ID>] -- the .species key is
    // how SET_DIFF_CALCULATIONS joins lineage_index against species_index. A
    // target_group_* glob emits a List<path> when >1 group matched, or a bare Path
    // when exactly 1 -- fan_out_target_group handles both.
    def fan_out_target_group = { meta, files ->
        def file_list = files instanceof List ? files : [files]
        file_list.collect { f -> [[ID: f.getParent().getName(), species: meta.ID], f] }
    }

    COLOR_MAPPING.out.target_group_file_colors
    | flatMap(fan_out_target_group)
    | set { lineage_file_colors }

    COLOR_MAPPING.out.target_group_label_mapping
    | flatMap(fan_out_target_group)
    | set { lineage_label_mapping }

    GGCAT_GROUP(lineage_file_colors)

    SBWT_BUILD_GROUP(GGCAT_GROUP.out.unitigs)

    SBWT_BUILD_GROUP.out.index
    | map { meta, sbwt, lcs -> [meta, sbwt] }
    | set { lineage_sbwt_only }

    SBWT_CHECK_GROUP(lineage_sbwt_only)

    SBWT_CHECK_GROUP.out.index
    | join(SBWT_BUILD_GROUP.out.index.map { meta, sbwt, lcs -> [meta, lcs] })
    | set { lineage_checked_index } // tuple(meta, sbwt, lcs), meta.species set -- Index B

    lineage_file_colors
    | join(lineage_checked_index)
    | set { lineage_themisto_build_input }

    THEMISTO_BUILD_GROUP(lineage_themisto_build_input)
    THEMISTO_STATS_GROUP(THEMISTO_BUILD_GROUP.out.index)
    THEMISTO_EXPORT_GROUP(THEMISTO_STATS_GROUP.out.index)


    // ===================================================================
    // Step 07 - candidate index filtering (per lineage): core_catchall_filter.py
    // filters that lineage's own export down to a candidate marker unitig FASTA,
    // using params.candidate_min_freq/candidate_min_genome_count. Rebuilt into a
    // real index below -- candidate_index (E) is Python-derived until then, not
    // yet a real SBWT index.
    // ===================================================================

    THEMISTO_EXPORT_GROUP.out.unitigs
    | join(THEMISTO_EXPORT_GROUP.out.color_sets)
    | join(THEMISTO_EXPORT_GROUP.out.metadata)
    | join(lineage_label_mapping)
    | set { candidate_filter_input } // tuple(meta, unitigs, color_sets, export_metadata, label_mapping)

    CANDIDATE_FILTER(candidate_filter_input)

    // If nothing cleared candidate filtering's thresholds, the FASTA is empty (0
    // bytes, per core_catchall_filter.py) -- skip rebuilding an index from nothing
    // for that lineage rather than handing GGCAT an empty input. log.warn (not just
    // the Python script's own stderr) so this is visible in the main pipeline log,
    // not just that task's own work-dir log file.
    CANDIDATE_FILTER.out.unitigs
    | filter { meta, fasta ->
        if (fasta.size() == 0) {
            log.warn("No candidate unitigs survived filtering for lineage '${meta.ID}' -- skipping candidate_index rebuild.")
            return false
        }
        true
    }
    // stage: 'candidate' -- meta.ID/species stay exactly the lineage's own (required
    // for SET_DIFF_CALCULATIONS' join), so this is what lets GGCAT/SBWT/THEMISTO_BUILD's
    // shared publishDir tell this rebuild's output apart from that lineage's own.
    | map { meta, fasta -> [meta + [stage: 'candidate'], fasta] }
    | set { candidate_fasta_nonempty }

    // Step 4a (per core_catchall_filter.py's docstring) - rebuild the filtered
    // candidate FASTA into a real index. CANDIDATE_COLOR_LIST wraps it as a one-line
    // colour-list file first -- GGCAT/THEMISTO_BUILD expect that shape (a list of
    // input file paths, one colour per line), not a raw sequences FASTA directly.
    CANDIDATE_COLOR_LIST(candidate_fasta_nonempty)

    GGCAT_CANDIDATE(CANDIDATE_COLOR_LIST.out.file_colors)
    SBWT_BUILD_CANDIDATE(GGCAT_CANDIDATE.out.unitigs)

    SBWT_BUILD_CANDIDATE.out.index
    | map { meta, sbwt, lcs -> [meta, sbwt] }
    | set { candidate_sbwt_only }

    SBWT_CHECK_CANDIDATE(candidate_sbwt_only)

    SBWT_CHECK_CANDIDATE.out.index
    | join(SBWT_BUILD_CANDIDATE.out.index.map { meta, sbwt, lcs -> [meta, lcs] })
    | set { candidate_checked_index } // tuple(meta, sbwt, lcs)

    // QC gate only, not a data source -- confirms Themisto2 can build a structurally
    // sound index from the rebuild and reports sane counts. No THEMISTO_EXPORT_CANDIDATE:
    // candidate_index only ever feeds sbwt difference downstream, which never touches
    // exported unitigs/color_sets, so exporting them would be wasted compute.
    CANDIDATE_COLOR_LIST.out.file_colors
    | join(candidate_checked_index)
    | set { candidate_themisto_build_input }

    THEMISTO_BUILD_CANDIDATE(candidate_themisto_build_input)
    THEMISTO_STATS_CANDIDATE(THEMISTO_BUILD_CANDIDATE.out.index)

    emit:
    // Only the SBWT indexes set_diff_calculations.nf actually needs. The species
    // Themisto2 build/export and COLOR_MAPPING's species label_mapping/stats still
    // run above (lineage_specificity_score.py, once wired in, needs the species
    // export) -- they're just not part of this subworkflow's public contract anymore.
    sbwt_index       = checked_index           // tuple(meta, sbwt, lcs) -- species-wide (A), feeds set_diff_calculations.nf
    lineage_index    = lineage_checked_index   // tuple(meta, sbwt, lcs), meta.species set -- feeds set_diff_calculations.nf as B
    candidate_index  = candidate_checked_index // tuple(meta, sbwt, lcs), meta.species set -- feeds set_diff_calculations.nf as E
}
