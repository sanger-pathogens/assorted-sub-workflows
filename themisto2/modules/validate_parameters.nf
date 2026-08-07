def validate_path_param(param_option, param, type = "file", mandatory = true) {
    def valid_types = ["file", "directory", "any"] // "any" = must exist, file or directory both fine
    if (!valid_types.any { it == type }) {
        throw new Exception("Invalid type '${type}'. Possibilities are ${valid_types}.")
    }

    def param_name = param_option - "--"
    if (param) {
        def file_param = file(param)
        if (!file_param.exists()) {
            throw new Exception("The given ${param_name} '${param}' does not exist.")
        } else if (
            (type == "file" && !file_param.isFile()) ||
            (type == "directory" && !file_param.isDirectory())
        ) {
            throw new Exception("The given ${param_name} '${param}' is not a ${type}.")
        }
    } else if (mandatory) {
        throw new Exception("No ${param_name} specified. Please specify one using the ${param_option} option.")
    }
}

def validate_parameters() {
    def errors = 0

    // --metadata is always required
    try {
        validate_path_param("--metadata", params.metadata)
    } catch (Exception e) {
        log.error(e.message)
        errors += 1
    }

    // --label_col is always required -- color_mapping.py needs a real column name, not a blank
    if (!params.label_col) {
        log.error("No --label_col specified. Please specify the metadata column to use as the label.")
        errors += 1
    }

    // --assembly_input is always required -- either a directory of assemblies or a
    // txt file listing assembly paths, one per line. Which kind it is gets sniffed
    // from the path itself later (in color_mapping.nf), not from which param was set.
    try {
        validate_path_param("--assembly_input", params.assembly_input, "any")
    } catch (Exception e) {
        log.error(e.message)
        errors += 1
    }

    if (errors > 0) {
        throw new Exception("${errors} error(s) detected in input parameters. See above for details.")
    }
}
