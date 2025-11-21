//
// Check input manifest and assign idenifier channels
// Expected to be used as submodule to input data to IRODS extractor
//

workflow IRODS_MANIFEST_PARSE {

    take:
    lane_manifest // file: /path/to/manifest_of_lanes.csv

    main:
    Channel
        .fromPath( lane_manifest )
        .ifEmpty {exit 1, "File is empty / Cannot find file at ${lane_manifest}"}
        .splitCsv ( header:true, strip:true, sep:',' )
        .map { createChannel(it) }
        .set { meta }

    emit:
    meta
}

def setMissingValue(String var, row, params) {
    if (row[var] == null) {
        // If manifest doesn't include column for given var, return CLI (params) var value
        params[var]
    } else {
        def value = row[var]?.toString()?.trim()
        if (value == "-1") {
            // If manifest includes column for given var, but set to -1, return -1
            -1
        } else {
            // Otherwise, if column left empty set to -1, else use value (string converted)
            value ?: -1
        }
    }
}

def createChannel(LinkedHashMap row) {
    def metadata = [:]
    
    // Convert fields to integers, defaulting to -1 if missing or empty
    ['studyid', 'runid', 'laneid', 'plexid'].each { key ->
        metadata[key] = row[key]?.toString()?.trim()?.isInteger() ? row[key].toInteger() : -1
    }
    
    // Ensure at least studyid or runid is provided if laneid or plexid is specified
    if (metadata.laneid != -1 || metadata.plexid != -1) {
        if (metadata.studyid == -1 && metadata.runid == -1) {
            log.warn ("Cannot submit an iRODS query based on laneid or plexid metadata tags where neither studyid or runid are specified, as this query would catch too many file objects.\nThe row ${row} in the input manifest is ignored")
            return "none"
        }
    }
    
    // Collect additional metadata fields
    def extraFields = row.keySet() - ['studyid', 'runid', 'laneid', 'plexid']
    extraFields.each { key ->
        if (row[key]?.toString()?.trim()) {
            metadata[key] = row[key].toString()
        }
    }
    
    return metadata
}
