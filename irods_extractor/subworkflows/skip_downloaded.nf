//Input: meta_cram_ch → tuple(meta, cram_path)
//Done check: scan params.outdir for FASTQ .gz using the same nested/flat glob; infer IDs from filenames; skip if metadata.ID is in inferred IDs
//Output: filtered channel same shape → tuple(meta, cram_path) (only those not already done)


workflow FILTER_EXISTING_OUTPUTS {
    take:
        meta_cram_ch // tuple(meta, cram_path)
    
    main:
        // Error check:
        if ( !params.outdir ) {
            error "FILTER_EXISTING_OUTPUTS: params.outdir must be set to use output-based skipping."
        } 

        // Build the glob once
        def fastq_glob = (params.save_method == "nested")
            ? "${params.outdir}/*/${params.preexisting_fastq_tag}/*${params.split_sep_for_ID_from_fastq}.gz"
            : "${params.outdir}/${params.preexisting_fastq_tag}/*${params.split_sep_for_ID_from_fastq}.gz"

        log.debug "FILTER_EXISTING_OUTPUTS: save_method=${params.save_method}; glob=${fastq_glob}"

        Channel.fromPath(fastq_glob).set{ preexisting_fastq_path_ch }

        // Count how many FASTQ files already exist at the start, and store that number so it can be reused later.
        preexisting_fastq_path_ch.count()
            .map { n ->
                log.info "FILTER_EXISTING_OUTPUTS: ${n} existing FASTQ files found in output directory."
                return n
        }
        .set { preexisting_fastq_count_ch }

       // Extract existing_id
        preexisting_fastq_path_ch.map{ preexisting_fastq_path ->
            preexisting_fastq_path.Name.split("${params.split_sep_for_ID_from_fastq}")[0]
        }
        .collect()
        .ifEmpty([])
        .set{ existing_id }
       
        // Extra debug + warning around ID inference quality
        existing_id.combine(preexisting_fastq_count_ch)
            .map { ids, n_fastqs ->
                log.debug "FILTER_EXISTING_OUTPUTS: ${ids.size()} existing sample IDs detected"

                if (n_fastqs > 0 && ids.size() == 0) {
                    log.warn "FILTER_EXISTING_OUTPUTS: Found ${n_fastqs} FASTQ files but inferred 0 IDs. Likely glob/ID-split mismatch (check split_sep_for_ID_from_fastq and output layout)."
                }

                // Optional debug: show first few inferred IDs (avoid spam)
                def preview_n = Math.min(5, ids.size())
                if (preview_n > 0) {
                    log.debug "FILTER_EXISTING_OUTPUTS: first ${preview_n} inferred IDs: ${ids.take(preview_n).join(', ')}"
                }

                return ids
            }
            .set { existing_id_checked }


        // Filter inputs against existing_id
        meta_cram_ch.combine( existing_id )
            .filter { metadata, cram_path, existing -> !(metadata.ID.toString() in existing) } 
            .map { it[0,1] }
            .set{ filtered_meta_cram_ch }
        // Debug count: how many will be processed:

        filtered_meta_cram_ch.toList().map { filtered_list ->
            def new_downloads = filtered_list.size()
            log.debug "FILTER_EXISTING_OUTPUTS: ${new_downloads} input items will be processed"
        }
        
    emit:
        filtered_meta_cram_ch

}