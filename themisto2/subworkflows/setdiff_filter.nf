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

workflow SET_DIFF_CALCULATIONS {
    take:
    /// tuple(meta, sbwt) -- e.g. ATB. meta.ID = "ATB"; broadcast to every species below.
    background_ch
    /// tuple(meta, sbwt, lcs) -- A: species/group-wide index. meta.ID = species id.
    /// From BUILD_COLOR_INDEX.out.sbwt_index (whole-species genome set).
    species_index_ch
    /// tuple(meta, sbwt, lcs) -- B: lineage-only index. meta.ID = lineage id,
    /// meta.species = parent species id (must match a species_index_ch meta.ID).
    /// PREREQUISITE, not built here -- see note at the bottom of this file.
    lineage_index_ch
    /// tuple(meta, sbwt, lcs) -- E: per-lineage candidate index (Step07-derived,
    /// built+checked upstream). meta.ID = lineage id, meta.species = parent species id.
    candidate_index_ch

    main:
    // bg_excl (step08 "C") = background - species. Species-agnostic, computed
    // once per species, reused by every lineage of that species below.
    background_ch
    | combine(species_index_ch)
    | map { bg_meta, bg_sbwt, sp_meta, sp_sbwt, sp_lcs ->
        [[ID: "bg_excl_${sp_meta.ID}", species: sp_meta.ID], bg_sbwt, sp_sbwt]
    }
    | set { bg_excl_input }

    SBWT_DIFFERENCE_BG_EXCL(bg_excl_input)
    SBWT_CHECK_BG_EXCL(SBWT_DIFFERENCE_BG_EXCL.out.index)

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
    SBWT_CHECK_BG_EXCL.out.index
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
    bg_excl  = SBWT_CHECK_BG_EXCL.out.index   // step08 "C" -- background_exclusion, see README glossary
    xlin_bg  = SBWT_CHECK_XLIN_BG.out.index   // step08 "D" -- cross_lineage_background
    lin_cand = SBWT_CHECK_LIN_CAND.out.index  // step08 "F" -- lineage_specific_candidates
    markers  = SBWT_CHECK_MARKERS.out.index   // step08 "G" -- final_candidate_markers
}

// TODO -- none of this subworkflow's take: channels have a real producer yet.
// main.nf includes BUILD_COLOR_INDEX but doesn't call SET_DIFF_CALCULATIONS at
// all. In call order:
//
// 1. background_ch (ATB/GTDB) -- not built by any subworkflow. Needs a param
//    for pre-built external .sbwt path(s), fed in via Channel.fromPath(...) in
//    main.nf. Doesn't exist yet.
//
// 2. lineage_index_ch (B) -- per 08_set_difference_filtering.md, get this by
//    calling BUILD_COLOR_INDEX again with metadata_ch/assembly_ch filtered to
//    one lineage's genomes (same mechanism as species-wide A, scoped smaller;
//    one real ggcat+sbwt build per lineage -- 765 GPSCs for s_pneumoniae, so
//    not free). Two blockers before that actually works:
//      - BUILD_COLOR_INDEX's emitted meta is just [ID: ...], no .species key,
//        but this channel needs meta.species to join() against species_index_ch.
//        Nothing injects that yet -- has to happen wherever this second call
//        is made.
//      - This is exactly what --target_groups/Index B (color_mapping.py) was
//        built for, but index_target_group/<group>/ isn't fed into
//        GGCAT/SBWT_BUILD anywhere -- BUILD_COLOR_INDEX still only ever builds
//        from index_species/, regardless of --target_groups.
//
// 3. candidate_index_ch (E) -- core_catchall_filter.py's output, rebuilt as a
//    real SBWT index (doc's "Step 4a"). No process wraps
//    core_catchall_filter.py as a Nextflow step yet, and no process rebuilds
//    its unitig-ID-list output via GGCAT+SBWT_BUILD.
//
// Also still out of scope, per the doc's own open question: G_gtdb
// (markers - a GTDB-based bg_excl) -- add a second background_ch/combine pair
// once GTDB is actually built and the "replaces vs. supplements ATB" question
// is settled, not before.
