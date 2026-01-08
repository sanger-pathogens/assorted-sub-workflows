// Input: meta_dataobj_ch -> tuple(meta, dataobj)
// Done check: scan params.outdir for FASTQ .gz using the same nested/flat glob;
// infer IDs from filenames; skip if metadata.ID is in inferred IDs
// Output: filtered channel same shape -> tuple(meta, dataobj) (only those not already done)

workflow FILTER_EXISTING_OUTPUTS {
    take:
        meta_dataobj_ch // tuple(meta, dataobj)

    main:
        // Build the glob once
        fastq_glob = (params.save_method == 'nested')
            ? "${params.outdir}/*/${params.preexisting_fastq_tag}/*${params.split_sep_for_ID_from_fastq}.gz"
            : "${params.outdir}/${params.preexisting_fastq_tag}/*${params.split_sep_for_ID_from_fastq}.gz"

        // delete: preexisting_fastq_path_ch = Channel.fromPath(fastq_glob, checkIfExists: true)
        preexisting_fastq_path_ch = Channel.fromPath(fastq_glob)


        /* Count how many FASTQ files already exist at the start
        preexisting_fastq_count_ch = preexisting_fastq_path_ch
            .count()
            .map { n ->
                log.info "FILTER_EXISTING_OUTPUTS: ${n} existing output files found in output directory."
                n
            }
        */
        
        // Extract IDs from filenames
        // collect() will emit once at completion; if there are no files we force it to emit []
        existing_id_ch = preexisting_fastq_path_ch
            .map { p -> p.getName().tokenize(params.split_sep_for_ID_from_fastq)[0] }
            .collect()
            .map { ids ->
                ids = (ids ?: []).unique()
                if (ids.isEmpty()) {
                    log.info "FILTER_EXISTING_OUTPUTS: No existing output file found - nothing will be skipped"
                }
                ids
            }
            .ifEmpty { [] } // <-- ensures channel emits once even if no files matched

        // Filter inputs against existing IDs
        do_not_exist_ch = meta_dataobj_ch
            .combine(existing_id_ch)
            .filter { meta, dataobj, existing_ids ->
                if (meta == null || meta.ID == null) {
                    log.DEBUG "FILTER_EXISTING_OUTPUTS: metadata.ID missing for item; NOT skipping. metadata=${meta}"
                    return true
                }
                !existing_ids.contains(meta.ID.toString())

            }
            .map { meta, dataobj, _ -> tuple(meta, dataobj) }
        /*
        do_not_exist_ch = meta_dataobj_ch
            .combine(existing_id_checked_ch)
            .filter { meta_cram, existing_ids ->
                def (meta, dataobj) = meta_cram
                if (meta == null || meta.ID == null) {
                    log.warn "FILTER_EXISTING_OUTPUTS: metadata.ID missing for item; NOT skipping. metadata=${meta}"
                    return true
                }
                !existing_ids.contains(meta.ID.toString())
            }
            .map { meta_cram, dataobj -> meta_cram }
        */
        // 🔎 DEBUG: check final structure
        do_not_exist_ch.view { meta, dataobj ->
            "FINAL → meta.ID=${meta?.ID}, dataobj=${dataobj}"
        }
        

    emit:
        do_not_exist_ch
}
