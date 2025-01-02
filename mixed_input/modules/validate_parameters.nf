def validate_path_param(
    param_option,
    param, 
    type = "file", 
    mandatory = true) {

    def valid_types = ["file", "directory"]
    if (!valid_types.any { it == type }) {
        throw new Exception("Invalid type '${type}'. Possibilities are ${valid_types}.")
        return 1
    }

    def param_name = (param_option - "--").replaceAll("", " ")
    if (param) {
        def file_param = file(param)
        if (!file_param.exists()) {
            throw new Exception("The given ${param_name} '${param}' does not exist.")
            return 1
        } else if (
            (type == "file" && !file_param.isFile()) ||
            (type == "directory" && !file_param.isDirectory())
        ) {
            throw new Exception("The given ${param_name} '${param}' is not a ${type}.")
            return 1
        }
    } else if (mandatory) {
        throw new Exception("No ${param_name} specified. Please specify one using the ${param_option} option.")
        return 1
    }
    return 0
}

def validate_integer(potentialInt) {
    /*
    Our param as a -1 is an int however when we give it e.g. 6495 on the command line it becomes a string
    this function takes a string or int (or any type) and attempts to turn it into a string then back to an int

    if it can it returns true if for any reason it cannot it returns false
    */
    if (potentialInt == null) return false
    try {
        Integer.parseInt(potentialInt.toString())
        return true
    } catch (NumberFormatException e) {
        return false
    }
}

def validate_parameters() {
    def errors = 0
    def workflows_to_run = []

    // Determine which input specification method is provided
    def manifest_ena_exists = params.manifest_ena != null
    def manifest_of_lanes_exists = params.manifest_of_lanes != null
    def manifest_of_reads_exists = params.manifest_of_reads != null
    def manifest_exists = params.manifest != null
    
    //and CLI
    def has_studyid = params.studyid != -1
    def has_runid = params.runid != -1
    def has_laneid = params.laneid != -1
    def has_plexid = params.plexid != -1

    // anything provided?
    def input_specified = manifest_ena_exists ||
                          manifest_of_lanes_exists ||
                          manifest_of_reads_exists ||
                          manifest_exists ||
                          has_studyid ||
                          has_runid ||
                          has_laneid ||
                          has_plexid

    if (!input_specified) {
        throw new Exception("""No input specification provided. Please specify one or a mix of:

                    Manifests:
                    - --manifest_ena
                    - --manifest_of_lanes
                    - --manifest_of_reads or --manifest

                    CLI Arguments:
                    - --studyid
                    - --runid
                    - --laneid
                    - --plexid""")
        errors += 1
    }

    // Validate manifest inputs
    if (manifest_of_lanes_exists) {
        errors += validate_path_param("--manifest_of_lanes", params.manifest_of_lanes)
        workflows_to_run << 'IRODS'
    }

    if (manifest_ena_exists) {
        errors += validate_path_param("--manifest_ena", params.manifest_ena)
        workflows_to_run << 'ENA'
    }

    if (manifest_exists) {
        if (manifest_of_reads_exists) {
            //exit early if you give both
            throw new Exception("Cannot supply both ${params.manifest_of_reads} and ${params.manifest} as they are alias's of the same manifest")
        } else {
            // If manifest is provided but manifest_of_reads is not, use manifest for reads
            log.info("No --manifest_of_reads provided. Using --manifest parameter value '${params.manifest}' as --manifest_of_reads.\n")
            /*
            The manifest is not replaced at at this stage in reality as that is handled in the mixed_input workflow
            def manifestToUse = params.manifest_of_reads ? params.manifest_of_reads : params.manifest
            This is just a good opportunity to throw a error
            */
        }

        errors += validate_path_param("--manifest", params.manifest)
        workflows_to_run << 'READS_MANIFEST'
    }

    if (manifest_of_reads_exists) {
        errors += validate_path_param("--manifest_of_reads", params.manifest_of_reads)
        workflows_to_run << 'READS_MANIFEST'
    }

    // Validate CLI-based inputs
    if (has_studyid || has_runid || has_laneid || has_plexid) {
        // Validate individual parameters
        if (has_studyid && !validate_integer(params.studyid)) {
            throw new Exception("Invalid studyid provided: ${params.studyid}")
            errors += 1
        }
        if (has_runid && !validate_integer(params.runid)) {
            throw new Exception("Invalid runid provided: ${params.runid}")
            errors += 1
        }
        if (has_laneid && !validate_integer(params.laneid)) {
            throw new Exception("Invalid laneid provided: ${params.laneid}")
            errors += 1
        }
        if (has_plexid && !validate_integer(params.plexid)) {
            throw new Exception("Invalid plexid provided: ${params.plexid}")
            errors += 1
        }
        
        workflows_to_run << 'IRODS'
    }

    if (errors > 0) {
        throw new Exception(String.format("%d errors detected", errors))
    }

    return workflows_to_run
}