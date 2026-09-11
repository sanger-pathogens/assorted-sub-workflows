// Everything downstream of the species-wide index build (build_color_index.nf):
//
//   1. LINEAGE-SPECIFICITY filtering -- lineage_specificity_filter.py picks candidate
//      marker unitigs out of the species-wide export (lineage-core, and lineage-specific
//      against sister lineages within the same species). Emits one candidate FASTA per
//      targeted lineage.
//   2. Candidate index rebuild -- candidate FASTA -> colour-list -> GGCAT -> SBWT ->
//      Themisto2 (QC/dump only, no export needed).
//   3. ATB CROSS-SPECIES check -- replaces the old bg_excl/markers sbwt set-diff (PAT-3570:
//      sbwt difference is colour-blind and doesn't scale at species-index level). Dumps the
//      candidate index to FASTA, pseudoaligns it against ATB-species.thm2
//      (THEMISTO2_ATB_PSEUDOALIGN, modules/themisto2.nf), then filter_atb_markers.py scores each candidate marker's
//      hit_frac against every ATB species colour and keeps only markers that are solidly
//      within the target species (>=atb_min_within) and essentially absent from every
//      other one (<=atb_max_outside). See filter_atb_markers.nf / bin/filter_atb_markers.py.
//
// A species with no atb_target_species (manifest column blank -- e.g. not present in ATB
// at all) skips step 3 for that species: its rebuilt candidate markers pass straight
// through UNCHECKED, with a loud log.warn, rather than being silently dropped or the whole
// run failing. Every other species goes through the full check.
include { LINEAGE_SPECIFICITY_FILTER; CANDIDATE_COLOR_LIST } from '../modules/lineage_specificity_filtering.nf'
include { GGCAT as GGCAT_CANDIDATE                         } from '../modules/ggcat.nf'
include { SBWT_BUILD as SBWT_BUILD_CANDIDATE; SBWT_CHECK as SBWT_CHECK_CANDIDATE; SBWT_DUMP_UNITIGS } from '../modules/sbwt.nf'
include { THEMISTO2_BUILD as THEMISTO2_BUILD_CANDIDATE; THEMISTO2_STATS as THEMISTO2_STATS_CANDIDATE; THEMISTO2_ATB_PSEUDOALIGN } from '../modules/themisto2.nf'
include { FILTER_ATB_MARKERS                               } from '../modules/filter_atb_markers.nf'
include { CHECKPOINT_COUNT                                 } from '../modules/checkpoint_count.nf'

workflow MARKER_FILTERING {
    take:
    /// tuple(meta, unitigs, color_sets, export_metadata, label_mapping) -- species-wide
    /// Themisto2 export, from BUILD_COLOR_INDEX.out.species_export. meta.ID = species id.
    species_export_ch
    /// tuple(meta, target_groups_string) -- MANIFEST_PARSE.out.target_groups
    target_groups_ch
    /// tuple(meta, atb_target_species_string) -- MANIFEST_PARSE.out.atb_target_species.
    /// Blank string = species skipped for the ATB check (see header comment above).
    atb_target_species_ch

    main:
    // Fixed ATB reference files -- staged once as value channels, broadcast to every task
    // by Nextflow (same pattern the old bg_index used, just plain path() inputs instead of
    // a channel built from Channel.fromPath since neither process needs meta on these).
    atb_index_ch       = Channel.value(file(params.atb_index, checkIfExists: true))
    atb_color_names_ch = Channel.value(file(params.atb_color_names, checkIfExists: true))

    // ============ 1. Lineage-specificity candidate filtering ============
    // Runs ONCE per species over the SPECIES-wide export. target_groups is joined in HERE
    // only, on the slim [ID: species] meta -- see manifest_parse.nf for why it's kept out
    // of meta upstream.
    species_export_ch
    | join(target_groups_ch)
    | set { specificity_filter_input } // tuple(meta, unitigs, color_sets, export_metadata, label_mapping, target_groups)

    LINEAGE_SPECIFICITY_FILTER(specificity_filter_input)

    // One task emits one candidate FASTA per requested lineage. Fan them into one item
    // each, meta = [ID: lineage, species: species run ID] -- the lineage ID comes from
    // the filename, meta.species is the join key for the rest of this subworkflow.
    LINEAGE_SPECIFICITY_FILTER.out.unitigs
    | flatMap { meta, files ->
        def file_list = files instanceof List ? files : [files]
        file_list.collect { f -> [[ID: (f.name - '_candidate_unitigs.fasta'), species: meta.ID], f] }
    }
    | set { candidate_fasta_per_lineage } // tuple(meta, fasta) -- one per lineage, may be empty

    // Empty FASTA = nothing cleared the thresholds -- skip the rebuild for that lineage.
    candidate_fasta_per_lineage
    | filter { meta, fasta ->
        if (fasta.size() == 0) {
            log.warn("No candidate unitigs survived filtering for lineage '${meta.ID}' -- skipping candidate_index rebuild.")
            return false
        }
        true
    }
    | map { meta, fasta -> [meta + [stage: 'candidate'], fasta] }
    | set { candidate_fasta_nonempty }

    // ============ 2. Candidate index rebuild ============
    // CANDIDATE_COLOR_LIST wraps the FASTA as a one-line colour-list file first --
    // GGCAT/THEMISTO2_BUILD expect that, not a raw FASTA.
    CANDIDATE_COLOR_LIST(candidate_fasta_nonempty)

    GGCAT_CANDIDATE(CANDIDATE_COLOR_LIST.out.file_colors)
    SBWT_BUILD_CANDIDATE(GGCAT_CANDIDATE.out.unitigs)

    SBWT_BUILD_CANDIDATE.out.index
    | map { meta, sbwt, lcs -> [meta, sbwt] }
    | set { candidate_sbwt_only }

    SBWT_CHECK_CANDIDATE(candidate_sbwt_only)

    SBWT_CHECK_CANDIDATE.out.index
    | join(SBWT_BUILD_CANDIDATE.out.index.map { meta, sbwt, lcs -> [meta, lcs] })
    | set { candidate_checked_index } // tuple(meta, sbwt, lcs), meta.species set

    // QC gate only -- confirms the rebuild is structurally sound. No export needed: the
    // ATB check below reads unitigs straight off SBWT_DUMP_UNITIGS, not a Themisto2 export.
    CANDIDATE_COLOR_LIST.out.file_colors
    | join(candidate_checked_index)
    | set { candidate_themisto_build_input }

    THEMISTO2_BUILD_CANDIDATE(candidate_themisto_build_input)
    THEMISTO2_STATS_CANDIDATE(THEMISTO2_BUILD_CANDIDATE.out.index)

    // Candidate index -> FASTA, needed as the ATB pseudoalign query (and as the fallback
    // "unchecked" output for species with no atb_target_species, below).
    candidate_checked_index
    | map { meta, sbwt, lcs -> [meta, sbwt] }
    | set { candidate_sbwt_for_dump }

    SBWT_DUMP_UNITIGS(candidate_sbwt_for_dump)

    // ============ 3. ATB cross-species check ============
    atb_target_species_ch
    | map { meta, atb -> [meta.ID, atb] } // meta.ID = species (this channel is species-level)
    | set { atb_target_species_by_species }

    SBWT_DUMP_UNITIGS.out.unitigs
    | map { meta, fasta -> [meta.species, meta, fasta] } // meta.species = join key (this channel is per-lineage)
    | join(atb_target_species_by_species)
    | map { species, meta, fasta, atb -> [meta, fasta, atb] }
    | branch {
        meta, fasta, atb ->
        checked: atb?.trim()
        unchecked: true
    }
    | set { atb_branch }

    atb_branch.unchecked
    | map { meta, fasta, atb ->
        log.warn("No atb_target_species set for species '${meta.species}' -- lineage '${meta.ID}' candidate markers pass through the ATB cross-species check UNVERIFIED.")
        [meta, fasta]
    }
    | set { markers_unchecked }

    atb_branch.checked
    | map { meta, fasta, atb -> [meta, fasta] }
    | set { atb_pseudoalign_input }

    THEMISTO2_ATB_PSEUDOALIGN(atb_pseudoalign_input, atb_index_ch)

    THEMISTO2_ATB_PSEUDOALIGN.out.jsonl
    | join(atb_pseudoalign_input)
    | join(atb_branch.checked.map { meta, fasta, atb -> [meta, atb] })
    | set { filter_atb_input } // tuple(meta, jsonl, fasta, atb_target_species)

    FILTER_ATB_MARKERS(filter_atb_input, atb_color_names_ch)

    // Final markers = PASS fasta for checked species, raw candidate fasta (unverified) for
    // species with no ATB mapping. Both keyed the same way: meta.ID = lineage, meta.species.
    FILTER_ATB_MARKERS.out.pass
    | mix(markers_unchecked)
    | set { markers }

    // ============ Checkpoint counts ============
    // `order` keys continue after build_color_index.nf's stages (which end at 40).
    Channel.empty()
    | mix( candidate_fasta_per_lineage.map           { meta, f -> [meta, 50, 'candidate_specificity_filter', 'fasta', f] } )
    | mix( THEMISTO2_BUILD_CANDIDATE.out.index.map   { meta, f -> [meta, 60, 'candidate_rebuilt_index', 'themisto', f] } )
    | mix( SBWT_DUMP_UNITIGS.out.unitigs.map         { meta, f -> [meta, 70, 'candidate_dumped_fasta', 'fasta', f] } )
    | mix( markers.map                               { meta, f -> [meta, 80, 'markers_atb_checked', 'fasta', f] } )
    | set { checkpoint_inputs }

    CHECKPOINT_COUNT(checkpoint_inputs)

    emit:
    markers     = markers                        // tuple(meta, fasta) -- final candidate markers, meta.species/meta.ID(lineage) set
    checkpoints = CHECKPOINT_COUNT.out.row        // tuple(meta, row_tsv) -- per-stage count rows
}
