// General helper functions
def validate_path_param(
    param_option, 
    param, 
    type="file", 
    mandatory=true) {
        valid_types=["file", "directory"]
        if (!valid_types.any { it == type }) {
                log.error("Invalid type '${type}'. Possibilities are ${valid_types}.")
                return 1
        }
        param_name = (param_option - "--").replaceAll("_", " ")
        if (param) {
            def file_param = file(param)
            if (!file_param.exists()) {
                log.error("The given ${param_name} '${param}' does not exist.")
                return 1
            } else if (
                (type == "file" && !file_param.isFile())
                ||
                (type == "directory" && !file_param.isDirectory())
            ) {
                log.error("The given ${param_name} '${param}' is not a ${type}.")
                return 1
            }
        } else if (mandatory) {
            log.error("No ${param_name} specified. Please specify one using the ${param_option} option.")
            return 1
        }
        return 0
    }

def validate_choice_param(param_option, param, choices) {
    param_name = (param_option - "--").replaceAll("_", " ")
    if (!choices.contains(param)) {
        log.error("Please specify the ${param_name} using the ${param_option} option. Possibilities are ${choices}.")
        return 1
    }
    return 0
}

def validate_number_param(param_option, param) {
    param_name = (param_option - "--").replaceAll("_", " ")
    if (param != null) /* Explicit comparison with null, because 0 is an acceptable value */ {
        if (!(param instanceof Number)) {
            log.error("The ${param_name} specified with the ${param_option} option must be a valid number")
            return 1
        }
    } else {
        log.error("Please specify the ${param_name} using the ${param_option} option")
        return 1
    }
    return 0
}

// Pipeline-specific helper functions
def validate_kraken2_db_param(param_option, param) {
    is_invalid_path = validate_path_param(param_option, param, type="directory", mandatory=true)
    if (is_invalid_path) {
        return 1
    }
    required_files = ["hash.k2d", "opts.k2d", "taxo.k2d"]
    kraken2_db_path = file(param)
    if (!required_files.any { file("${kraken2_db_path}/${it}").isFile() }) {
        log.error("The Kraken2 database provided to ${param_option} does not contain one of the required files: ${required_files}.")
        return 1
    }
    return 0
}

def validate_kraken2_db_for_bracken_build(param_option, param) {
    validate_kraken2_db_param(param_option, param)
    kraken2_db_dir = file(param, checkIfExists: true)
    bracken_kmer_distrib = file("${kraken2_db_dir}/database${params.read_len}mers.kmer_distrib")
    if (bracken_kmer_distrib.exists()) {
        return 0
    }
    required_dirs = ["library"]
    if (!required_dirs.any { file("${kraken2_db_path}/${it}").isDirectory() }) {
        log.error("The Kraken2 database provided to ${param_option} does not contain one of the required directories for building a bracken database: ${required_dirs}.")
        return 1
    }
    return 0
}

def validate_parameters() {
    def errors = 0

    errors += validate_path_param("--input", params.input)
    errors += validate_kraken2_db_param("--kraken2_db", params.kraken2_db)
    errors += validate_kraken2_db_for_bracken_build("--kraken2_db", params.kraken2_db)
    errors += validate_number_param("--read_len", params.read_len)
    errors += validate_number_param("--kmer_len", params.kmer_len)
    errors += validate_choice_param("--classification_level", params.classification_level, ['D','P','C','O','F','G','S'])
    errors += validate_number_param("--threshold", params.threshold)
    errors += validate_choice_param("--get_classified_reads", params.get_classified_reads, [true, false])
    errors += validate_number_param("--kraken2_threads", params.kraken2_threads)
    errors += validate_number_param("--bracken_threads", params.bracken_threads)
    errors += validate_choice_param("--enable_building", params.enable_building, [true, false])
    
    if (errors > 0) {
        log.error(String.format("%d errors detected", errors))
        exit 1
    }
}
