include { LINEAGE_SPECIFICITY_FILTER; CANDIDATE_COLOUR_LIST } from '../modules/lineage_specificity_filtering.nf'
include { GGCAT as GGCAT_CANDIDATE                         } from '../modules/ggcat.nf'
include { SBWT_BUILD as SBWT_BUILD_CANDIDATE; SBWT_CHECK as SBWT_CHECK_CANDIDATE; SBWT_DUMP_UNITIGS as SBWT_DUMP_UNITIGS_CANDIDATE } from '../modules/sbwt.nf'
include { THEMISTO2_BUILD as THEMISTO2_BUILD_CANDIDATE; THEMISTO2_STATS as THEMISTO2_STATS_CANDIDATE; THEMISTO2_EXPORT as THEMISTO2_EXPORT_CANDIDATE; THEMISTO2_ATB_PSEUDOALIGN } from '../modules/themisto2.nf'
include { FILTER_ATB_MARKERS                               } from '../modules/filter_atb_markers.nf'
include { CHECKPOINT_FASTA; CHECKPOINT_THEMISTO             } from '../modules/checkpoint.nf'

workflow MARKER_FILTERING {
    take:
    species_export_ch
    target_groups_ch
    atb_target_species_ch
    atb_exclude_species_ch

    main:
    atb_index_ch       = Channel.value(file(params.atb_index, checkIfExists: true))
    atb_colour_names_ch = Channel.value(file(params.atb_colour_names, checkIfExists: true))

    species_export_ch
    | join(target_groups_ch)
    | set { specificity_filter_input }

    LINEAGE_SPECIFICITY_FILTER(specificity_filter_input)

    LINEAGE_SPECIFICITY_FILTER.out.unitigs
    | flatMap { meta, files ->
        def file_list = files instanceof List ? files : [files]
        file_list.collect { f -> [[ID: (f.name - '_candidate_unitigs.fasta'), species: meta.ID], f] }
    }
    | set { candidate_fasta_per_lineage }

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

    CANDIDATE_COLOUR_LIST(candidate_fasta_nonempty)

    GGCAT_CANDIDATE(CANDIDATE_COLOUR_LIST.out.file_colours)
    SBWT_BUILD_CANDIDATE(GGCAT_CANDIDATE.out.unitigs)

    SBWT_BUILD_CANDIDATE.out.index
    | map { meta, sbwt, lcs -> [meta, sbwt] }
    | set { candidate_sbwt_only }

    SBWT_CHECK_CANDIDATE(candidate_sbwt_only)

    // Pass or fail goes to .nextflow.log; a failed index goes no further (QC-only rebuild,
    // so this only loses the candidate_sbwt/themisto/export checkpoint rows for the group).
    SBWT_CHECK_CANDIDATE.out.index
    | filter { meta, sbwt, result ->
        if (result == 'PASS') {
            log.info("SBWT check passed: ${meta.species} candidate index ${meta.ID}")
            return true
        }
        log.warn("SBWT check FAILED: ${meta.species} candidate index ${meta.ID} -- skipping its Themisto2 rebuild and SBWT dump.")
        false
    }
    | map { meta, sbwt, result -> [meta, sbwt] }
    | join(SBWT_BUILD_CANDIDATE.out.index.map { meta, sbwt, lcs -> [meta, lcs] })
    | set { candidate_checked_index }

    CANDIDATE_COLOUR_LIST.out.file_colours
    | join(candidate_checked_index)
    | set { candidate_themisto_build_input }

    THEMISTO2_BUILD_CANDIDATE(candidate_themisto_build_input)
    THEMISTO2_STATS_CANDIDATE(THEMISTO2_BUILD_CANDIDATE.out.index)

    // Checkpoint only: one record per unitig (no strand duplicates), unlike SBWT_DUMP_UNITIGS_CANDIDATE.
    THEMISTO2_EXPORT_CANDIDATE(THEMISTO2_STATS_CANDIDATE.out.index)

    candidate_checked_index
    | map { meta, sbwt, lcs -> [meta, sbwt] }
    | set { candidate_sbwt_for_dump }

    SBWT_DUMP_UNITIGS_CANDIDATE(candidate_sbwt_for_dump)

    atb_target_species_ch
    | join(atb_exclude_species_ch)
    | map { meta, atb, excl -> [meta.ID, atb, excl] }
    | set { atb_target_species_by_species }

    // ATB check reads GGCAT's candidate unitigs directly: one record per unitig. The candidate
    // SBWT/Themisto2 rebuild is QC-only (checkpoint stages candidate_sbwt_unitigs,
    // candidate_themisto_index, candidate_export_unitigs); SBWT_DUMP_UNITIGS_CANDIDATE writes
    // both strands of every unitig, which used to double the ATB query and its counts.
    GGCAT_CANDIDATE.out.unitigs
    | map { meta, fasta -> [meta.species, meta, fasta] }
    | join(atb_target_species_by_species)
    | map { species, meta, fasta, atb, excl -> [meta, fasta, atb, excl] }
    | branch {
        meta, fasta, atb, excl ->
        checked: atb?.trim()
        unchecked: true
    }
    | set { atb_branch }

    atb_branch.unchecked
    | map { meta, fasta, atb, excl ->
        log.warn("Species '${meta.species}' isn't in ATB -- lineage '${meta.ID}' candidate markers pass through the ATB cross-species check UNVERIFIED.")
        [meta, fasta]
    }
    | set { markers_unchecked }

    atb_branch.checked
    | map { meta, fasta, atb, excl -> [meta, fasta] }
    | set { atb_pseudoalign_input }

    THEMISTO2_ATB_PSEUDOALIGN(atb_pseudoalign_input, atb_index_ch)

    THEMISTO2_ATB_PSEUDOALIGN.out.jsonl
    | join(atb_pseudoalign_input)
    | join(atb_branch.checked.map { meta, fasta, atb, excl -> [meta, atb, excl] })
    | set { filter_atb_input }

    FILTER_ATB_MARKERS(filter_atb_input, atb_colour_names_ch)

    FILTER_ATB_MARKERS.out.pass
    | mix(markers_unchecked)
    | set { markers }

    // Stage names and their order in the report live in checkpoint.nf (checkpoint_steps()).
    Channel.empty()
    | mix( candidate_fasta_per_lineage.map              { meta, f -> [meta, 'candidate_specificity_filter', 'fasta', f] } )
    | mix( GGCAT_CANDIDATE.out.unitigs.map              { meta, f -> [meta, 'candidate_ggcat_unitigs', 'fasta', f] } )
    | mix( SBWT_DUMP_UNITIGS_CANDIDATE.out.unitigs.map  { meta, f -> [meta, 'candidate_sbwt_unitigs', 'fasta', f] } )
    | mix( THEMISTO2_EXPORT_CANDIDATE.out.unitigs.map   { meta, f -> [meta, 'candidate_export_unitigs', 'fasta', f] } )
    | mix( markers.map                                  { meta, f -> [meta, 'markers_atb_checked', 'fasta', f] } )
    | set { fasta_checkpoint_inputs }

    THEMISTO2_STATS_CANDIDATE.out.stats.map { meta, f -> [meta, 'candidate_themisto_index', f] }
    | set { themisto_checkpoint_inputs }

    CHECKPOINT_FASTA(fasta_checkpoint_inputs)
    CHECKPOINT_THEMISTO(themisto_checkpoint_inputs)

    emit:
    markers     = markers
    checkpoints = CHECKPOINT_FASTA.out.row.mix(CHECKPOINT_THEMISTO.out.row)
}
