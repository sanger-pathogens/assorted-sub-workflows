def validate_path_param(
    paramoption, 
    param, 
    type = "file", 
    mandatory = true) {

    def valid_types = ["file", "directory"]
    if (!valid_types.any { it == type }) {
        log.error("Invalid type '${type}'. Possibilities are ${valid_types}.")
        return 1
    }

    def param_name = (paramoption - "--").replaceAll("", " ")
    if (param) {
        def file_param = file(param)
        if (!file_param.exists()) {
            log.error("The given ${param_name} '${param}' does not exist.")
            return 1
        } else if (
            (type == "file" && !file_param.isFile()) ||
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

def validate_parameters() {
    def errors = 0
    def workflows_to_run = []

    // Determine which input specification method is provided
    def manifest_of_lanes_exists = params.manifest_of_lanes != null
    def manifest_of_reads_exists = params.manifest_of_reads != null
    def manifest_exists = params.manifest != null
    
    //and CLI
    def has_studyid = params.studyid != -1
    def has_runid = params.runid != -1
    def has_laneid = params.laneid != -1
    def has_plexid = params.plexid != -1

    // anything provided?
    def input_specified = manifest_of_lanes_exists || 
                          manifest_of_reads_exists || 
                          manifest_exists || 
                          has_studyid || 
                          has_runid || 
                          has_laneid || 
                          has_plexid

    if (!input_specified) {
    log.error("""No input specification provided. Please specify one of:

                Manifests:
                - --manifest_of_lanes
                - --manifest_of_reads or --manifest

                CLI Identifiers:
                - --studyid
                - --runid
                - --laneid
                - --plexid""")
    errors += 1
} 

    // Validate and route based on input type
    if (manifest_of_lanes_exists || manifest_of_reads_exists || manifest_exists) {
        // Validate manifests
        if (manifest_of_lanes_exists) {
            errors += validate_path_param("--manifest_of_lanes", params.manifest_of_lanes)
            workflows_to_run << 'IRODS'
        }

        if (manifest_exists) {
            // If manifest is provided but manifest_of_reads is not, use manifest for reads
            if (!manifest_of_reads_exists) {
                log.info("manifest_of_reads not provided. Using manifest as manifest_of_reads.")
                params.manifest_of_reads = params.manifest
            }
            
            errors += validate_path_param("--manifest", params.manifest)
            workflows_to_run << 'READS_MANIFEST'
        }

        if (manifest_of_reads_exists) {
            errors += validate_path_param("--manifest_of_reads", params.manifest_of_reads)
            workflows_to_run << 'READS_MANIFEST'
        }
    }

    // Check for CLI-based inputs
    if (has_studyid || has_runid || has_laneid || has_plexid) {
        // Validate individual parameters
        if (has_studyid && params.studyid == -1) {
            log.error("Invalid studyid provided.")
            errors += 1
        }
        if (has_runid && params.runid == -1) {
            log.error("Invalid runid provided.")
            errors += 1
        }
        if (has_laneid && params.laneid == -1) {
            log.error("Invalid laneid provided.")
            errors += 1
        }
        if (has_plexid && params.plexid == -1) {
            log.error("Invalid plexid provided.")
            errors += 1
        }
        
        workflows_to_run << 'IRODS'
    }

    if (errors > 0) {
        log.error(String.format("%d errors detected", errors))
        System.exit(1)
    }

    return workflows_to_run
}