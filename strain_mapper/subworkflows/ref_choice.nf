include { REF_MANIFEST_PARSE } from './ref_manifest.nf'

workflow CHECK_REFERENCES {
    main:
    generic_reference = null
    reference_manifest = null

    if (params.reference) {
        generic_reference = file(params.reference, checkIfExists: true)
    } else {
        if (!params.reference_manifest) {
            log.error "You must provide either a reference fasta file (option `--reference`) or a reference manifest (option `--reference_manifest`)."
            exit 1
        } else {
            log.info "No generic reference provided (option `--reference`), will use reference manifest to determine references for each sample. Samples without a reference in the manifest will be skipped."
        }
    }

    if (params.reference_manifest) {
        reference_manifest = file(params.reference_manifest, checkIfExists: true)
    }

    ch_reference_manifest = channel.empty()

    if (reference_manifest) {
        REF_MANIFEST_PARSE(reference_manifest)

        REF_MANIFEST_PARSE.out.references
            .map { metaref, reference -> [metaref.ID, metaref, reference] }
            .set { ch_reference_manifest }
    }

    emit:
    generic_reference
    ch_reference_manifest
}

workflow PICK_REFERENCE {
    take:
    reads_ch // tuple( meta, read_1, read_2 )
    ch_reference_manifest // tuple( ID, meta, reference path )
    generic_reference // value: Path

    main:
    reads_ch
    .map { metaread, reads_1, reads_2 ->
        [metaread.ID, metaread, reads_1, reads_2]
    }
    .join(ch_reference_manifest, remainder: true)
    .map { row ->
        def (mid, meta, reads_1, reads_2) = row
        def reference = row.size() >= 6 ? row[5] : null
        [meta, reads_1, reads_2, reference ?: generic_reference] // prefer manifest reference if available
    }
    .set { ch_reads_with_ref }

    if (!params.drop_without_ref){
        ch_reads_with_ref.filter { meta, reads_1, reads_2, reference -> reference == null }
        .map { meta, reads_1, reads_2, reference -> meta.ID }
        .toList()
        .subscribe { ids ->
            if (ids) {
                log.error "No reference could be assigned to ${ids.size()} sample(s): ${ids.sort().join(', ')}\n" +
                        "Add them to --reference_manifest, supply a fallback with --reference, " +
                        "or set --drop_without_ref to skip them instead."
                exit 1
            }
        }            
    }

    ch_reads_with_ref.filter { meta, reads_1, reads_2, reference ->
        // meta is null for reference manifest rows whose ID matched no input sample,
        // reference is null for samples with no manifest entry and no --reference
        meta != null && reference != null
    }
    .set { all_reads_ready_to_map_with_ref_ch }

    emit:
    all_reads_ready_to_map_with_ref_ch
}