include { GUBBINS; GUBBINS_MASK } from './modules/gubbins.nf' 
include { SNP_SITES } from './modules/snp-sites.nf'
include { RAXML_NG } from './modules/raxml-ng.nf'


workflow CONSTRUCT_PHYLO {
    take:
    input_msa
    
    main:
    // Mask regions where recombination is likely to have occurred
    if (params.remove_recombination) {
        GUBBINS(input_msa)
        | GUBBINS_MASK

        GUBBINS_MASK.out.masked_msa.set { msa }
    } else {
        input_msa.set { msa }
    }

    // Infer subsitution model
    MODEL_FINDER(msa)

    // Extract SNP sites, get constant site frequencies and build tree
    SNP_SITES(msa)
    | RAXML_NG

    emit:
    BUILD_TREE.out.tree_channel
}