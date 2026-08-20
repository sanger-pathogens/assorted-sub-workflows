// SBWT_DIFFERENCE/SBWT_CHECK are generic (diff two SBWT indexes, verify the
// result) -- aliased once per stage here so the same two building blocks run
// every step08 stage without a separate process definition per stage. See
// 08_set_difference_filtering.md for the full methodology this mirrors, and
// the README's glossary for what each short name below stands for.
//
//   stage    step08  formula
//   bg_excl  C       background - species_index (A)
//   xlin_bg  D       species_index (A) - lineage_index (B)
//   lin_cand F       candidate_index (E) - xlin_bg
//   markers  G       lin_cand - bg_excl
include { SBWT_DIFFERENCE as SBWT_DIFFERENCE_BG_EXCL; SBWT_CHECK as SBWT_CHECK_BG_EXCL } from '../modules/sbwt.nf'
include { SBWT_DIFFERENCE as SBWT_DIFFERENCE_XLIN_BG; SBWT_CHECK as SBWT_CHECK_XLIN_BG } from '../modules/sbwt.nf'
include { SBWT_DIFFERENCE as SBWT_DIFFERENCE_LIN_CAND; SBWT_CHECK as SBWT_CHECK_LIN_CAND } from '../modules/sbwt.nf'
include { SBWT_DIFFERENCE as SBWT_DIFFERENCE_MARKERS; SBWT_CHECK as SBWT_CHECK_MARKERS } from '../modules/sbwt.nf'
// Not a step08 diff stage -- just verifies --bg_excl_index (an arbitrary,
// unverified user-supplied file used directly AS bg_excl -- see the bg_excl
// block below) actually loads before it's used downstream.
include { SBWT_CHECK as SBWT_CHECK_BACKGROUND } from '../modules/sbwt.nf'

workflow SET_DIFF_CALCULATIONS {
    take:
    // (no background_ch here -- it's built from params.bg_excl_index / params.bg_index
    // in main: below, not threaded through by the caller. See that block's comment.)
    /// tuple(meta, sbwt, lcs) -- A: species/group-wide index. meta.ID = species id.
    /// From BUILD_COLOR_INDEX.out.sbwt_index (whole-species genome set).
    species_index_ch
    /// tuple(meta, sbwt, lcs) -- B: lineage-only index. meta.ID = lineage id,
    /// meta.species = parent species id (must match a species_index_ch meta.ID).
    /// From BUILD_COLOR_INDEX.out.lineage_index -- only non-empty when
    /// --target_groups is set (see build_color_index.nf).
    lineage_index_ch
    /// tuple(meta, sbwt, lcs) -- E: per-lineage candidate index (Step07-derived,
    /// built+checked upstream). meta.ID = lineage id, meta.species = parent species id.
    candidate_index_ch

    main:
    // bg_excl (step08 "C") = background - species. --bg_excl_index lets the user
    // reuse an already-built bg_excl directly and skip re-running the hugemem-scale
    // diff against --bg_index (the raw background, e.g. ATB); if not given, it's
    // computed for real below. Either way this converges on bg_excl_checked, one
    // entry per species -- in the --bg_index branch, combine(species_index_ch) is
    // what makes the diff PER SPECIES (not one global diff); the --bg_excl_index
    // branch just tags the reused file with each species_index_ch entry's ID the
    // same way, since there's no diff left to run.
    if (params.bg_excl_index) {
        // bg_excl_index is an arbitrary user-supplied file, never verified by this
        // pipeline -- check it loads before using it downstream. Tag with the
        // file's own basename (e.g. ATB_all) instead of a generic label, so
        // publishDir paths/output filenames/logs identify which background
        // file was actually used.
        def bg_excl_name = file(params.bg_excl_index).baseName

        SBWT_CHECK_BACKGROUND(
            Channel.fromPath(params.bg_excl_index, checkIfExists: true)
                .map { sbwt -> [[ID: "bg_excl_${bg_excl_name}"], sbwt] }
        )

        SBWT_CHECK_BACKGROUND.out.index
        | combine(species_index_ch)
        | map { bg_meta, bg_sbwt, sp_meta, sp_sbwt, sp_lcs ->
            [[ID: "bg_excl_${bg_excl_name}_${sp_meta.ID}", species: sp_meta.ID], bg_sbwt]
        }
        | set { bg_excl_checked }
    } else {
        Channel.fromPath(params.bg_index, checkIfExists: true)
        | map { sbwt -> [[ID: 'bg_index'], sbwt] }
        | combine(species_index_ch)
        | map { bg_meta, bg_sbwt, sp_meta, sp_sbwt, sp_lcs ->
            [[ID: "bg_excl_${sp_meta.ID}", species: sp_meta.ID], bg_sbwt, sp_sbwt]
        }
        | set { bg_excl_input }

        SBWT_DIFFERENCE_BG_EXCL(bg_excl_input)
        SBWT_CHECK_BG_EXCL(SBWT_DIFFERENCE_BG_EXCL.out.index)
        SBWT_CHECK_BG_EXCL.out.index | set { bg_excl_checked }
    }

    // xlin_bg (step08 "D") = species-wide index (A) - lineage-only index (B)
    species_index_ch
    | map { meta, sbwt, lcs -> [meta.ID, sbwt] }
    | set { species_by_id }

    lineage_index_ch
    | map { meta, sbwt, lcs -> [meta.species, meta, sbwt] }
    | join(species_by_id)
    | map { species, lin_meta, lin_sbwt, sp_sbwt ->
        [[ID: "xlin_bg_${lin_meta.ID}", species: species, lineage: lin_meta.ID], sp_sbwt, lin_sbwt]
    }
    | set { xlin_bg_input }

    SBWT_DIFFERENCE_XLIN_BG(xlin_bg_input)
    SBWT_CHECK_XLIN_BG(SBWT_DIFFERENCE_XLIN_BG.out.index)

    // lin_cand (step08 "F") = candidate index (E) - xlin_bg
    SBWT_CHECK_XLIN_BG.out.index
    | map { meta, sbwt -> [meta.lineage, sbwt] }
    | set { xlin_bg_by_lineage }

    candidate_index_ch
    | map { meta, sbwt, lcs -> [meta.ID, meta, sbwt] }
    | join(xlin_bg_by_lineage)
    | map { lineage, cand_meta, cand_sbwt, d_sbwt ->
        [[ID: "lin_cand_${lineage}", species: cand_meta.species, lineage: lineage], cand_sbwt, d_sbwt]
    }
    | set { lin_cand_input }

    SBWT_DIFFERENCE_LIN_CAND(lin_cand_input)
    SBWT_CHECK_LIN_CAND(SBWT_DIFFERENCE_LIN_CAND.out.index)

    // markers (step08 "G") = lin_cand - bg_excl
    bg_excl_checked
    | map { meta, sbwt -> [meta.species, sbwt] }
    | set { bg_excl_by_species }

    SBWT_CHECK_LIN_CAND.out.index
    | map { meta, sbwt -> [meta.species, meta, sbwt] }
    | join(bg_excl_by_species)
    | map { species, f_meta, f_sbwt, c_sbwt ->
        [[ID: "markers_${f_meta.lineage}", species: species, lineage: f_meta.lineage], f_sbwt, c_sbwt]
    }
    | set { markers_input }

    SBWT_DIFFERENCE_MARKERS(markers_input)
    SBWT_CHECK_MARKERS(SBWT_DIFFERENCE_MARKERS.out.index)

    emit:
    bg_excl  = bg_excl_checked                // step08 "C" -- background_exclusion, see README glossary
    xlin_bg  = SBWT_CHECK_XLIN_BG.out.index   // step08 "D" -- cross_lineage_background
    lin_cand = SBWT_CHECK_LIN_CAND.out.index  // step08 "F" -- lineage_specific_candidates
    markers  = SBWT_CHECK_MARKERS.out.index   // step08 "G" -- final_candidate_markers
}

// STATUS: all inputs this subworkflow needs now have real producers -- the three
// take: channels (species_index_ch/lineage_index_ch/candidate_index_ch) come from
// BUILD_COLOR_INDEX.out.sbwt_index/.lineage_index/.candidate_index (wired in
// main.nf; meta.species is set correctly on both B and E for the join()s above),
// and bg_excl (C) is built from params.bg_excl_index/params.bg_index inside this
// subworkflow itself -- neither is a take: input. lineage_index_ch/
// candidate_index_ch (B/E) are only non-empty when --target_groups is set.
//
// Still out of scope, per 08_set_difference_filtering.md's own open question:
// G_gtdb (markers - a GTDB-based bg_excl) -- add a second background/combine
// pair once GTDB is actually built and the "replaces vs. supplements ATB"
// question is settled, not before.
