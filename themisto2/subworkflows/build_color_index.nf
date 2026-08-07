include { COLOR_MAPPING                   } from '../modules/color_mapping.nf'
include { GGCAT                           } from '../modules/ggcat.nf'
include { SBWT_BUILD; SBWT_CHECK          } from '../modules/sbwt.nf'
include { THEMISTO_BUILD; THEMISTO_STATS; THEMISTO_EXPORT } from '../modules/themisto2.nf'
include { validate_parameters             } from '../modules/validate_parameters.nf'

workflow BUILD_COLOR_INDEX {
    take:
    /// Input files: metadata and assemblies
    metadata_ch
    assembly_ch

    main:
    validate_parameters()


    // TODO: only works for a single run -- combine() cross-multiplies if there's ever
    // >1 metadata/assembly. Switch to a manifest channel of pre-paired tuples for multi-run.
    
    // Tag metadata with an ID and pair with assembly input (plumbing only -- COLOR_MAPPING
    // needs one tuple, not two separate inputs).
    metadata_ch
    | map { metadata -> [["ID": metadata.baseName], metadata] }
    | combine(assembly_ch)
    | set { color_mapping_input }

    // Step 02 - map sample metadata + assemblies to Themisto's colour-file format
    COLOR_MAPPING(color_mapping_input)


    // Step 03 - build unitigs from the colour file
    GGCAT(COLOR_MAPPING.out.file_colors)

    // Step 04 - build the SBWT index from the unitigs, then verify it loads correctly
    // (verification split into its own process -- much lighter workload than the build,
    // no reason to hold a cpu_16/mem_32 reservation just to run a quick check)
    SBWT_BUILD(GGCAT.out.unitigs)
    SBWT_CHECK(SBWT_BUILD.out.index)


    // Step 05 - build the Themisto2 index (needs the colour file *and* the verified SBWT index)
    COLOR_MAPPING.out.file_colors
    | join(SBWT_CHECK.out.index)
    | set { themisto_build_input }

    THEMISTO_BUILD(themisto_build_input)
    THEMISTO_STATS(THEMISTO_BUILD.out.index)

    // Step 06 - export the index
    THEMISTO_EXPORT(THEMISTO_STATS.out.index)

    emit:
    index            = THEMISTO_STATS.out.index      // tuple(meta, index_thm2)
    export_unitigs   = THEMISTO_EXPORT.out.unitigs
    color_sets       = THEMISTO_EXPORT.out.color_sets
    export_metadata  = THEMISTO_EXPORT.out.metadata
    label_mapping    = COLOR_MAPPING.out.label_mapping
    stats            = COLOR_MAPPING.out.stats
}
