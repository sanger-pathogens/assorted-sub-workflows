include { MINIBWA_INDEX;
          MINIBWA                 } from '../assemble/modules/minibwa.nf'
include { INDEX;
          SORT_TO_BAM             } from './modules/samtools.nf'
include { CONTIG_DEPTHS_NO_INTRA  } from './modules/metabat2.nf'
include { SPLIT_DEPTHS;
          MAXBIN2                 } from './modules/maxbin2.nf'
include { COMEBIN                 } from './modules/comebin.nf'
include { SEMIBIN2                } from './modules/semibin2.nf'
include { METACAT                 } from './modules/metacat.nf'

/*
##############################################################################################################################################################
#
# NOTE: as of this fork, METABAT1/METABAT2/CONCOCT and their modules (metabat2.nf's METABAT1/METABAT2, concoct.nf)
# are no longer wired into MAG_BINNING below - superseded by comebin + SemiBin2 + MetaCAT + MaxBin2, a combo
# validated (18-sample benchmark, exhaustive 3-/4-binner combo sweep) to statistically match the full original
# 6-binner ensemble (Wilcoxon p=0.41) while being ~44% cheaper end-to-end, and to beat the original metaWRAP-trio-only
# pipeline by +22.8% total HQ bins (p=0.0006). The module files themselves are left untouched/importable in case
# a future comparison wants them back.
#
##############################################################################################################################################################
*/

/*
##############################################################################################################################################################
#
# This script is meant to be run on the outputs of assembly.sh pipeline to split the assembly contigs into metagenomic bins.
# Ideally it should take in the assembly file of all of your samples, followed by the reads of all the samples that went into the assembly.
# The more samples, the better the binning. 
#
# The script uses metaBAT2 and/or CONCOCT and/or MaxBin2 to bin the contigs. MetaBAT2 is the default due to its speed and great performance,
# but all these binners have their advantages and disadvantages, so it recommended to run the bin_refinement module to QC the bins, get the 
# best bins of all of each method, and to reassembly and refine the final bins. 
#
# Author of pipeline: German Uritskiy. I do not clain any authorship of the many programs this pipeline uses.
# For questions, bugs, and suggestions, contact me at guritsk1@jhu.edu.
#
# Modified by Sam Dougan into nextflow :)
##############################################################################################################################################################
*/

// Function to join multiple binning channels into a single one
// that emits: [meta, [bin1, bin2, bin3]]
def join_bins(List binning_channels) {
    // Start with the first channel
    def joined_channel = binning_channels[0]

    // Iteratively join with the remaining channels
    (1..<binning_channels.size()).each { index ->
        def next_channel = binning_channels[index]

        joined_channel = joined_channel.join(next_channel)
            .map { tuple1, tuple2 ->
                def meta = tuple1[0]  // shared metadata (e.g. sample ID)

                // Get list of current bins from previous step
                def bin_list = (tuple1[1] instanceof List) ? tuple1[1] : [tuple1[1]]

                // Add the new bin to the list
                def new_bin = tuple2[1]
                return [meta, bin_list + [new_bin]]
            }
    }

    return joined_channel
}

workflow MAG_BINNING {
    take:
    contigs
    reads

    main:

    MINIBWA_INDEX(contigs)
    | set { indexed_contigs }

    reads.join(indexed_contigs)
    | MINIBWA
    | SORT_TO_BAM
    | set { bam }

    INDEX(bam)
    | set { bam_plus_index }

    bam_plus_index
    | join(contigs)
    | set { bam_bai_and_contigs }

    COMEBIN_WF(bam_bai_and_contigs)
    | set { comebin_bins }

    SEMIBIN2_WF(bam_bai_and_contigs)
    | set { semibin2_bins }

    METACAT_WF(bam_bai_and_contigs)
    | set { metacat_bins }

    MAXBIN_WF(bam, contigs)
    | set { maxbin_bins }

    comebin_bins
    | join(semibin2_bins)
    | join(metacat_bins)
    | join(maxbin_bins)
    | set { final_bins }

    emit:
    final_bins
}

workflow COMEBIN_WF {
    take:
    bam_bai_and_contigs  // [meta, bam, bai, assembly]

    main:
    COMEBIN(bam_bai_and_contigs)
    | set { bins }

    emit:
    bins
}

workflow SEMIBIN2_WF {
    take:
    bam_bai_and_contigs  // [meta, bam, bai, assembly]

    main:
    SEMIBIN2(bam_bai_and_contigs)
    | set { bins }

    emit:
    bins
}

workflow METACAT_WF {
    take:
    bam_bai_and_contigs  // [meta, bam, bai, assembly]

    main:
    METACAT(bam_bai_and_contigs)
    | set { bins }

    emit:
    bins
}

workflow MAXBIN_WF {
    take:
    bam
    contigs

    main:
    CONTIG_DEPTHS_NO_INTRA(bam)
    | SPLIT_DEPTHS
    | join(contigs)
    | MAXBIN2
    | set { bins }

    emit:
    bins
}
