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

def createChannel(LinkedHashMap row) {
    def metadata = [:]
    
    // Convert fields to integers, defaulting to -1 if missing or empty
    ['studyid', 'runid', 'laneid', 'plexid'].each { key ->
        metadata[key] = row[key]?.toString()?.trim()?.isInteger() ? row[key].toInteger() : -1
    }
    
    // Ensure at least studyid or runid is provided if laneid or plexid is specified
    if (metadata.laneid != -1 || metadata.plexid != -1) {
        if (metadata.studyid == -1 && metadata.runid == -1) {
            log.warn "Cannot submit an iRODS query with only laneid or plexid without studyid or runid. Row ${row} is ignored."
            return "none"
        }
    }
    
    // Collect additional metadata fields
    def extraFields = row.keySet() - ['studyid', 'runid', 'laneid', 'plexid', 'target', 'type']
    extraFields.each { key ->
        if (row[key]?.toString()?.trim()) {
            metadata[key] = row[key].toString()
        }
    }
    
    // Set target and type if necessary
    if (metadata.studyid != -1 || metadata.runid != -1 || extraFields.any { row[it]?.toString()?.trim() }) {
        metadata.target = row.target?.toString()?.trim() ?: "${params.target}"
        metadata.type = row.type?.toString()?.trim() ?: "${params.type}"
    } else if (row.target?.toString()?.trim() || row.type?.toString()?.trim()) {
        log.warn "Cannot submit an iRODS query with only target or type. Row ${row} is ignored."
        return "none"
    }
    
    return metadata
}
