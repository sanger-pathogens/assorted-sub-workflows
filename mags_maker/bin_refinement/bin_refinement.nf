include { BINETTE } from './modules/binette.nf'

/*
##############################################################################################################################################################
#
# Reconciles the four binners' outputs (comebin, SemiBin2, MetaCAT, MaxBin2) with Binette into one final bin set per sample.
#
# This replaces the original pipeline's custom n-wise bin-combination + BINNING_REFINER + CheckM1/2 + MERGE_BINS logic
# (see git history for that implementation) with a single, purpose-built bin-reconciliation tool. Binette does the
# same job - find the best-scoring combination of contigs across all input binners' candidate bins - more directly
# and was benchmarked (18-sample CHAIN+BBS cohort) to beat the original metaWRAP-trio-only pipeline reconciled the
# same way by +22.8% total HQ bins (Wilcoxon p=0.0006), and to statistically match reconciling all 6 original binners
# at once while costing ~44% less total pipeline time (comebin+SemiBin2+MetaCAT+MaxBin2 vs all 6, p=0.41).
#
##############################################################################################################################################################
*/

workflow MAG_BIN_REFINEMENT {
    take:
    bins      // [meta, comebin_bins, semibin2_bins, metacat_bins, maxbin2_bins]
    contigs   // [meta, assembly]

    main:
    bins
    | join(contigs)
    | BINETTE
    | set { binette_out }

    emit:
    final_bins = binette_out.results
}
