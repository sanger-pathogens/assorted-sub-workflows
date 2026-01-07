//Input: meta_cram_ch → tuple(meta, cram_path)
//Done check: scan params.outdir for FASTQ .gz using the same nested/flat glob; infer IDs from filenames; skip if metadata.ID is in inferred IDs
//Output: filtered channel same shape → tuple(meta, cram_path) (only those not already done)

workflow FILTER_EXISTING_OUTPUTS {
    take:
        meta_cram_ch // tuple(meta, cram_path)

    main:
        // Build the glob once
        def fastq_glob = (params.save_method == "nested")
            ? "${params.outdir}/*/${params.preexisting_fastq_tag}/*${params.split_sep_for_ID_from_fastq}.gz"
            : "${params.outdir}/${params.preexisting_fastq_tag}/*${params.split_sep_for_ID_from_fastq}.gz"
            
        //log.debug "FILTER_EXISTING_OUTPUTS: save_method=${params.save_method}; glob=${fastq_glob}"

        // Use normal Groovy variable assignment (avoids scope collisions)
        def preexisting_fastq_path_ch = Channel.fromPath(fastq_glob)

        /*
        preexisting_fastq_path_ch
            .take(5)
            .view { p -> "DEBUG preexisting_fastq_path_ch item: ${p} (class=${p?.getClass()?.name})" }
        */

        // Count how many FASTQ files already exist at the start
        def preexisting_fastq_count_ch = preexisting_fastq_path_ch
            .count()
            .map { n ->
                log.info "FILTER_EXISTING_OUTPUTS: ${n} existing FASTQ files found in output directory."
                n
            }
            
        /*
        preexisting_fastq_count_ch
            .view { n -> "DEBUG preexisting_fastq_count_ch: ${n} (class=${n?.getClass()?.name})" }
        */
        
        // Extract IDs from filenames (define once, then branch)
        def extracted_id_ch = preexisting_fastq_path_ch
            .map { p -> p.getName().split("${params.split_sep_for_ID_from_fastq}")[0] }

        /* debug branch (safe)
        extracted_id_ch
            .take(10)
            .view { id -> "DEBUG extracted_id (pre-collect): ${id} (class=${id?.getClass()?.name})" }
        */
        
        // real branch
        def existing_id_ch = extracted_id_ch
            .collect()
            .map { (it ?: []).unique() }
            

        // Extra debug + warning around ID inference quality
        def existing_id_checked_ch = existing_id_ch
            .combine(preexisting_fastq_count_ch)
            .map { ids, n_fastqs ->
                log.info "FILTER_EXISTING_OUTPUTS: ${n_fastqs} existing FASTQ files found"
                log.info "FILTER_EXISTING_OUTPUTS: ${ids.size()} existing sample IDs detected"

                if (n_fastqs > 0 && ids.size() == 0) {
                    log.warn "FILTER_EXISTING_OUTPUTS: Found ${n_fastqs} FASTQ files but inferred 0 IDs -> glob/ID split mismatch"
                }

                def preview_n = Math.min(5, ids.size())
                if (preview_n > 0) {
                    log.debug "FILTER_EXISTING_OUTPUTS: first ${preview_n} inferred IDs: ${ids.take(preview_n).join(', ')}"
                }

                ids
            }

        // Filter inputs against existing IDs
        def do_not_exist = meta_cram_ch
            .combine(existing_id_checked_ch)
            .filter { meta_cram, existing_ids ->
                def metadata = meta_cram[0]
                def id = metadata?.ID?.toString()

                if ( !id ) {
                    log.warn "FILTER_EXISTING_OUTPUTS: metadata.ID missing for item; NOT skipping. metadata=${metadata}"
                    return true
                }
                
                return !(id in existing_ids)
            }
            .map { meta_cram, _ -> meta_cram }


        log.warn "FILTER_EXISTING_OUTPUTS LOADED: emitting do_not_exist"

    emit:
        do_not_exist
}
