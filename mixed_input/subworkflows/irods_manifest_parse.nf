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
        .map { create_channel(it) }
        .set { meta }

    emit:
    meta
}

def create_channel(LinkedHashMap row) {
    def meta = [:]
    meta.studyid = ((! row.studyid) || ("${row.studyid}" == "")) ? -1 : "${row.studyid}".toInteger()
    meta.runid = ((! row.runid) || ("${row.runid}" == "")) ? -1 : "${row.runid}".toInteger()
    meta.laneid = ((! row.laneid) || ("${row.laneid}" == "")) ? -1 : "${row.laneid}".toInteger()
    meta.plexid = ((! row.plexid) || ("${row.plexid}" == "")) ? -1 : "${row.plexid}".toInteger()
    if (((meta.studyid == -1) && (meta.runid == -1)) && ((meta.laneid != -1) || (meta.plexid != -1))) {
        log.warn ("Cannot submit an iRODS query based on laneid or plexid metadata tags where neither studyid or runid are specified, as this query would catch too many file objects.\nThe row ${row} in the input manifest is ignored")
        return "none"
    }
    def extraFields = row.keySet().minus(['studyid', 'runid', 'laneid', 'plexid', 'target', 'type'])
    extraFields.each { key ->
        if ("${row[key]}" != "") {
            meta[key] = row[key].toString()
        }
    }
    if ((meta.studyid != -1) || (meta.runid != -1) || (extraFields.any { "${row[it]}" != "" })) {
        meta.target = ((! row.target) || ("${row.target}" == "")) ? "1" : "${row.target}"
        meta.type = ((! row.type) || ("${row.type}" == "")) ? "cram" : "${row.type}"
    }
    else if ((row.target && ("${row.target}" != "")) || (row.type && (row.type != ""))) {
        log.warn ("Cannot submit an iRODS query solely based on target or type metadata tags, as this query would catch too many file objects.\nThe row ${row} in the input manifest is ignored")
        return "none"
    }
    return meta
}
