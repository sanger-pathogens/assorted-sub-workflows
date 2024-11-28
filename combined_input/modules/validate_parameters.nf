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
    log.error("""No input specification provided. Please specify one or a mix of:

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

    // Validate manifests if provided
    if (manifest_of_lanes_exists) {
        errors += validate_path_param("--manifest_of_lanes", params.manifest_of_lanes)
    }

    if (manifest_exists) {
        // If manifest is provided but manifest_of_reads is not, use manifest for reads
        if (!manifest_of_reads_exists) {
            log.info("manifest_of_reads not provided. Using manifest as manifest_of_reads.")
            /*
            You can't actually overwrite the param.manifest_of_lanes at this stage to prevent having
            to return and pick this up later in the COMBINE_READS stage we do this
            def manifestToUse = params.manifest_of_reads ? params.manifest_of_reads : params.manifest
            */
        }
        
        errors += validate_path_param("--manifest", params.manifest)
    }

    if (manifest_of_reads_exists) {
        errors += validate_path_param("--manifest_of_reads", params.manifest_of_reads)
    }

    //todo move validation for the CLI methods from the workflow - clear errors if giving plex but nothing else for example

    if (errors > 0) {
        log.error(String.format("%d errors detected please correct and run again", errors))
        System.exit(1)
    }

    emit:
    Channel.of("Input parameters validated")
}