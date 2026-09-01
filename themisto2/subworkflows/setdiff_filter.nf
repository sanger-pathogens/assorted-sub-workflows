// SBWT_DIFFERENCE/SBWT_CHECK are generic (diff two SBWT indexes, verify the
// result) -- aliased once per step08 stage. See 08_set_difference_filtering.md
// and the README glossary for the full methodology.
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
// Verifies a user-supplied --bg_excl_index loads before it's used below.
include { SBWT_CHECK as SBWT_CHECK_BACKGROUND } from '../modules/sbwt.nf'
include { CHECKPOINT_COUNT_SBWT } from '../modules/checkpoint_count.nf'

workflow SET_DIFF_CALCULATIONS {
    take:
    // bg_excl (C) isn't a take: input -- built below from params.bg_excl_index/bg_index.
    /// tuple(meta, sbwt, lcs) -- A: species/group-wide index. meta.ID = species id.
    species_index_ch
    /// tuple(meta, sbwt, lcs) -- B: lineage-only index. meta.ID = lineage id,
    /// meta.species = parent species id. Only non-empty when --target_groups is set.
    lineage_index_ch
    /// tuple(meta, sbwt, lcs) -- E: per-lineage candidate index. meta.ID = lineage id,
    /// meta.species = parent species id.
    candidate_index_ch

    main:
    // bg_excl (C) = background - species. --bg_excl_index reuses a pre-built one and
    // skips the hugemem --bg_index diff. Both branches -> bg_excl_checked, one per species.
    if (params.bg_excl_index) {
        // Tag with the file basename so output names show which background was used.
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

    // xlin_bg (D) = species-wide index (A) - lineage-only index (B)
    species_index_ch
    | map { meta, sbwt, lcs -> [meta.ID, sbwt] }
    | set { species_by_id }

    lineage_index_ch
    | map { meta, sbwt, lcs -> [meta.species, meta, sbwt] }
    | join(species_by_id)
    | map { species, lin_meta, lin_sbwt, sp_sbwt ->
        // ID carries species + lineage so the diff_index filename identifies both
        [[ID: "xlin_bg_${species}_${lin_meta.ID}", species: species, lineage: lin_meta.ID], sp_sbwt, lin_sbwt]
    }
    | set { xlin_bg_input }

    SBWT_DIFFERENCE_XLIN_BG(xlin_bg_input)
    SBWT_CHECK_XLIN_BG(SBWT_DIFFERENCE_XLIN_BG.out.index)

    // lin_cand (F) = candidate index (E) - xlin_bg
    SBWT_CHECK_XLIN_BG.out.index
    | map { meta, sbwt -> [meta.lineage, sbwt] }
    | set { xlin_bg_by_lineage }

    candidate_index_ch
    | map { meta, sbwt, lcs -> [meta.ID, meta, sbwt] }
    | join(xlin_bg_by_lineage)
    | map { lineage, cand_meta, cand_sbwt, d_sbwt ->
        [[ID: "lin_cand_${cand_meta.species}_${lineage}", species: cand_meta.species, lineage: lineage], cand_sbwt, d_sbwt]
    }
    | set { lin_cand_input }

    SBWT_DIFFERENCE_LIN_CAND(lin_cand_input)
    SBWT_CHECK_LIN_CAND(SBWT_DIFFERENCE_LIN_CAND.out.index)

    // markers (G) = lin_cand - bg_excl
    bg_excl_checked
    | map { meta, sbwt -> [meta.species, sbwt] }
    | set { bg_excl_by_species }

    SBWT_CHECK_LIN_CAND.out.index
    | map { meta, sbwt -> [meta.species, meta, sbwt] }
    | join(bg_excl_by_species)
    | map { species, f_meta, f_sbwt, c_sbwt ->
        [[ID: "markers_${species}_${f_meta.lineage}", species: species, lineage: f_meta.lineage], f_sbwt, c_sbwt]
    }
    | set { markers_input }

    SBWT_DIFFERENCE_MARKERS(markers_input)
    SBWT_CHECK_MARKERS(SBWT_DIFFERENCE_MARKERS.out.index)

    // ============ Checkpoint counts ============
    // Distinct k-mers surviving each set-diff stage (see build_color_index.nf).
    // Side-channel only. stage keys continue the funnel order after A/B/E.
    Channel.empty()
    | mix( bg_excl_checked.map               { meta, sbwt -> [meta, 'C_bg_excl_08_diff', sbwt] } )
    | mix( SBWT_CHECK_XLIN_BG.out.index.map  { meta, sbwt -> [meta, 'D_xlin_bg_08_diff', sbwt] } )
    | mix( SBWT_CHECK_LIN_CAND.out.index.map { meta, sbwt -> [meta, 'F_lin_cand_08_diff', sbwt] } )
    | mix( SBWT_CHECK_MARKERS.out.index.map  { meta, sbwt -> [meta, 'G_markers_08_diff', sbwt] } )
    | set { checkpoint_inputs }

    CHECKPOINT_COUNT_SBWT(checkpoint_inputs)

    emit:
    bg_excl  = bg_excl_checked                // C -- background_exclusion
    xlin_bg  = SBWT_CHECK_XLIN_BG.out.index   // D -- cross_lineage_background
    lin_cand = SBWT_CHECK_LIN_CAND.out.index  // F -- lineage_specific_candidates
    markers  = SBWT_CHECK_MARKERS.out.index   // G -- final_candidate_markers
    checkpoints = CHECKPOINT_COUNT_SBWT.out.row // tuple(meta, row_tsv) -- per-stage count rows
}

// TODO: G_gtdb (markers - a GTDB-based bg_excl) -- add once GTDB is built.
