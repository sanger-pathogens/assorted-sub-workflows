//Input: meta_cram_ch → tuple(meta, cram_path)
//Done check: scan opts.outdir for FASTQ .gz (or caller-supplied outdir) ...
//Output: filtered channel same shape → tuple(meta, cram_path) (only those not already done)


workflow FILTER_EXISTING_OUTPUTS {
    take:
        meta_cram_ch // tuple(meta, cram_path)
        opts //opts is just a Groovy Map passed at call time
    
    main:
        opts = opts ?: [:]

        def outdir  = opts.outdir
        def tag     = opts.preexisting_fastq_tag ?: 'fastqs'
        def sep     = opts.split_sep_for_ID_from_fastq ?: '_1.fastq'
        def method  = opts.save_method ?: 'flat'

        if( !outdir ) error "FILTER_EXISTING_OUTPUTS: opts.outdir must be provided"

        log.debug "FILTER_EXISTING_OUTPUTS: outdir=${outdir} tag=${tag} sep=${sep} method=${method}"

        // Build the glob once
        def fastq_glob = (method == "nested")
            ? "${outdir}/*/${tag}/*${sep}.gz"
            : "${outdir}/${tag}/*${sep}.gz"

        log.debug "FILTER_EXISTING_OUTPUTS: save_method=${method}; glob=${fastq_glob}"



        def preexisting_fastq_path_ch = Channel.fromPath(fastq_glob)

// Count + ID extraction in one place (avoid `.count()` Integer stream)
// Reason: `.count()` emits a single Integer (e.g. 0) which can be mis-handled downstream,
// especially when combined with other channels. Collecting IDs gives us both:
// - ids.size() = FASTQ count
// - ids.unique() = list of inferred sample IDs

        def existing_ids_ch =
                    preexisting_fastq_path_ch
                        .map { p -> p.getName().split("${sep}")[0] }
                        .collect()
                        .map { ids ->
                            // ids is a List<String> of inferred IDs (may be empty)
                            log.info "FILTER_EXISTING_OUTPUTS: ${ids.size()} existing FASTQ files found in output directory."

                            def uniq = ids.unique()
                            log.debug "FILTER_EXISTING_OUTPUTS: ${uniq.size()} existing sample IDs detected"

                            if( ids.size() > 0 && uniq.size() == 0 ) {
                                log.warn "FILTER_EXISTING_OUTPUTS: Found ${ids.size()} FASTQ files but inferred 0 IDs. Check opts.split_sep_for_ID_from_fastq and output layout."
                            }

                            def preview_n = Math.min(5, uniq.size())
                            if( preview_n > 0 ) {
                                log.debug "FILTER_EXISTING_OUTPUTS: first ${preview_n} inferred IDs: ${uniq.take(preview_n).join(', ')}"
                            }

                            uniq
                        }
                        .ifEmpty { [] }   // ALWAYS a List, even when no FASTQs exist

// ---- Filter inputs ---------------------------------------------------
// combine shape: (tuple(meta, cram_path), existing_ids_list)
// Unpack inside closure to avoid arg-shape issues in DSL2.
        def filtered_meta_cram_ch =
            meta_cram_ch
                .combine(existing_ids_ch)
                .filter { meta_cram, existing_ids ->
                    def (metadata, cram_path) = meta_cram
                    def id = metadata.ID?.toString()
                    id && !(id in existing_ids)
                }
                .map { meta_cram, existing_ids -> meta_cram }

// ---- Optional lightweight debug (does not materialise whole stream) --
// If you need more logging, uncomment:
// filtered_meta_cram_ch.view { "FILTER_EXISTING_OUTPUTS: keeping ${it[0]?.ID} -> ${it[1]}" }

    emit:
        filtered_meta_cram_ch
}