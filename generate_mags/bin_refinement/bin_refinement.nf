include { BINNING_REFINER               } from './modules/bin_refiner.nf'
include { CHECKM2;
          CHECKM2 as CHECKM2_MERGED_BIN;
          CHECKM2 as CHECKM2_FINAL_BIN  } from './modules/checkm2.nf'
include { CHECKM;
          SUMMARISE_CHECKM              } from './modules/checkm1.nf'
include { MERGE_BINS                    } from './modules/merge_bins.nf'
include { DEREPLICATE_CONTIGS           } from './modules/dereplicate_contigs.nf'

/*

##############################################################################################################################################################
#
# This script is meant to be run on the outputs of binning.sh pipeline to analyze the metagenomic bins and arrive at the best possible putative genomes.
# 
# Author of pipeline: German Uritskiy. I do not clain any authorship of the many programs this pipeline uses.
# For questions, bugs, and suggestions, contact me at guritsk1@jhu.edu.
# 
# Modified by Yan Shao 1) checkM threshold default to 50/5 to align with medium-quality MAGs criteria; 2) clean-up bin folders  
##############################################################################################################################################################
*/

// Recursive function to get all combinations of size `n` from a list `fullBinList`
def getCombinations(List fullBinList, int n) {
    //exit early if small
    if (n == 0) return [[]]
    if (fullBinList.isEmpty() || n > fullBinList.size()) return []

    def firstBin = fullBinList[0]

    // Remaining fullBinList after the first (or empty list)
    def remainingBins = fullBinList.size() > 1 ? fullBinList[1..-1] : []

    // Include the first item in the combination
    def combosWithFirst = getCombinations(remainingBins, n - 1).collect { combo ->
        [firstBin] + combo
    }

    // Exclude the first item from the combination
    def combosWithoutFirst = getCombinations(remainingBins, n)

    // All valid combinations are either with or without the first item
    return combosWithFirst + combosWithoutFirst
}


workflow METAWRAP_BIN_REFINEMENT {
    take:
    bins
    reads

    main:

    bins
    | map { whole_channel -> 
        def meta = whole_channel[0]
        def bin_list = whole_channel[1..-1] //everything after meta
        return [meta, bin_list]
    }
    | set { all_bins }

    all_bins
    | transpose
    | map { meta, bin -> 
        def bin_name = bin.name
        [ meta, bin_name, bin ]
    }
    | set{ individual_bins }

    all_bins
    | flatMap { meta, bin_list ->
        // Get all 2-wise, 3-wise, ..., n-wise combinations
        def combinations = (2..bin_list.size()).collectMany { n ->
            getCombinations(bin_list, n)
        }

        // Map to new meta & bin group
        combinations.collect { combo ->
            def bin_names = combo.collect { path -> path.name } 
            def bins_label = bin_names.join('-')
            return [ meta, bins_label, combo ]
        }
    }
    | set { bin_combos }

    BINNING_REFINER(bin_combos)
    | mix(individual_bins)
    | set { refined_bins }
    
    if (params.checkm1) {
        CHECKM(refined_bins)
        | SUMMARISE_CHECKM
        | set { checkm }
    } else {
        CHECKM2(refined_bins)

        CHECKM2.out.results
        | set { checkm }
    }
    
    checkm
    | groupTuple
    | MERGE_BINS

    CHECKM2_MERGED_BIN(MERGE_BINS.out.merged_bins)
    | DEREPLICATE_CONTIGS
    | CHECKM2_FINAL_BIN

    CHECKM2_FINAL_BIN.out.results
    | set { final_bins }
    
    emit:
    final_bins
}