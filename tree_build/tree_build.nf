include { GUBBINS; GUBBINS_MASK } from './modules/gubbins.nf' 
include { SNP_SITES } from './modules/snp-sites.nf'
include { RAXML_NG } from './modules/raxml-ng.nf'
include { MODEL_FINDER } from './modules/iqtree.nf'
include { PLOT_TREE } from '../shared/modules/plotting.nf'

workflow CONSTRUCT_PHYLO {
    take:
    input_msa
    
    main:
    input_msa
    | map { msa_path -> ["ID": msa_path.baseName] }
    | set { meta }

    // Mask regions where recombination is likely to have occurred
    if (params.remove_recombination) {
        GUBBINS(input_msa)
        | GUBBINS_MASK

        GUBBINS_MASK.out.masked_msa.set { msa }
    } else {
        input_msa.set { msa }
    }

    // Infer substitution model
    if (params.check_model) {
        MODEL_FINDER(msa)
    }

    // Extract SNP sites, get constant site frequencies and build tree
    SNP_SITES(msa)
    | RAXML_NG

    meta
    | combine(RAXML_NG.out.tree)
    | set { plot_tree_input }

    PLOT_TREE(plot_tree_input)

    emit:
    RAXML_NG.out.tree
}