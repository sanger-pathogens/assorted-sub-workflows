// SBWT_DIFFERENCE/SBWT_CHECK are generic (diff two SBWT indexes, verify the
// result) -- aliased once per set-difference stage. See the stage table in the
// README for the full methodology.
//
//   stage    formula
//   bg_excl  background - species_index
//   markers  candidate_index - bg_excl
//
// The old cross-lineage set-diffs (xlin_bg = species_index - lineage_index;
// lin_cand = candidate_index - xlin_bg) were removed in PAT-3570: sbwt difference
// is colour-blind, so subtracting a lineage index from the species index never
// removes a k-mer a lineage shares with a sister lineage. Cross-lineage
// specificity is now a differential-frequency filter (lineage_specificity_filter.py),
// upstream of the candidate index.
include { SBWT_DIFFERENCE as SBWT_DIFFERENCE_BG_EXCL; SBWT_CHECK as SBWT_CHECK_BG_EXCL } from '../modules/sbwt.nf'
include { SBWT_DIFFERENCE as SBWT_DIFFERENCE_MARKERS; SBWT_CHECK as SBWT_CHECK_MARKERS } from '../modules/sbwt.nf'
// Verifies a user-supplied --bg_excl_index loads before it's used below.
include { SBWT_CHECK as SBWT_CHECK_BACKGROUND } from '../modules/sbwt.nf'
include { CHECKPOINT_COUNT_SBWT } from '../modules/checkpoint_count.nf'

workflow SET_DIFF_CALCULATIONS {
    take:
    // bg_excl isn't a take: input -- built below from params.bg_excl_index/bg_index.
    /// tuple(meta, sbwt, lcs) -- species/group-wide index. meta.ID = species id.
    species_index_ch
    /// tuple(meta, sbwt, lcs) -- per-lineage candidate index. meta.ID = lineage id,
    /// meta.species = parent species id.
    candidate_index_ch

    main:
    // bg_excl = background - species. --bg_excl_index reuses a pre-built one and
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

    // markers = candidate index - bg_excl, joined on the parent species
    bg_excl_checked
    | map { meta, sbwt -> [meta.species, sbwt] }
    | set { bg_excl_by_species }

    candidate_index_ch
    | map { meta, sbwt, lcs -> [meta.species, meta, sbwt] }
    | join(bg_excl_by_species)
    | map { species, cand_meta, cand_sbwt, c_sbwt ->
        [[ID: "markers_${species}_${cand_meta.ID}", species: species, lineage: cand_meta.ID], cand_sbwt, c_sbwt]
    }
    | set { markers_input }

    SBWT_DIFFERENCE_MARKERS(markers_input)
    SBWT_CHECK_MARKERS(SBWT_DIFFERENCE_MARKERS.out.index)

    // ============ Checkpoint counts ============
    // Distinct k-mers surviving each set-diff stage (see build_color_index.nf).
    // Side-channel only. `order` keys continue after the build stages (which end at 60).
    Channel.empty()
    | mix( bg_excl_checked.map               { meta, sbwt -> [meta, 70, 'bg_excl_diff', sbwt] } )
    | mix( SBWT_CHECK_MARKERS.out.index.map  { meta, sbwt -> [meta, 80, 'markers_diff', sbwt] } )
    | set { checkpoint_inputs }

    CHECKPOINT_COUNT_SBWT(checkpoint_inputs)

    emit:
    bg_excl  = bg_excl_checked                // background_exclusion
    markers  = SBWT_CHECK_MARKERS.out.index   // final_candidate_markers
    checkpoints = CHECKPOINT_COUNT_SBWT.out.row // tuple(meta, row_tsv) -- per-stage count rows
}

// TODO: a GTDB-based markers stage (markers - a GTDB-based bg_excl) -- add once GTDB is built.
