include { GUBBINS; GUBBINS_MASK } from '../modules/gubbins.nf' 
include { SNP_SITES } from '../modules/snp-sites.nf' 
include { BUILD_TREE } from '../modules/raxml-ng.nf'


workflow CONSTRUCT_PHYLO {
    take:
    input_msa        // msa
    
    main:
    if (params.remove_recombination) {
        GUBBINS(input_msa)
        | GUBBINS_MASK
        
        GUBBINS_MASK.out.masked_msa.set{ msa }
    } else {
        input_msa.set{ msa }
    }
    SNP_SITES(msa)
    | BUILD_TREE

    emit:
        BUILD_TREE.out.tree_channel
}