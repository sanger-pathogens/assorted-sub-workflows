def validate_parameters() {
    def errors = 0

    // --group_label / --sample_col are global (same metadata column names across the
    // whole run) and always required -- color_mapping.py needs real column names, not
    // blanks. The per-species metadata file and assembly path come in on the samples
    // channel now (see manifest_parse.nf), so they're validated by the caller, not here.
    if (!params.group_label) {
        log.error("No --group_label specified. Please specify the metadata column to use as the group label.")
        errors += 1
    }
    if (!params.sample_col) {
        log.error("No --sample_col specified. Please specify the metadata column holding the sample identifiers.")
        errors += 1
    }

    if (errors > 0) {
        throw new Exception("${errors} error(s) detected in input parameters. See above for details.")
    }
}
