/*
========================================================================================
    IMPORT MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULES
//

include { SKETCH_ASSEMBLIES; SKETCH_CORE_ACC_DIST } from './modules/sketchlib.nf'
include { RAPIDNJ                                 } from './modules/rapidnj.nf'
include { PLOT_TREE                               } from './modules/plotting.nf'


workflow SKETCHLIB_TREE {
    take:
    assemblies // [[ID: samples][assembly1, assembly2, ...]] collected channel of all assemblies to make tree from and a meta to name it

    main:
    SKETCH_ASSEMBLIES(assemblies)
    | SKETCH_CORE_ACC_DIST
    | RAPIDNJ
    | PLOT_TREE

    emit:
    RAPIDNJ.out.tree
}