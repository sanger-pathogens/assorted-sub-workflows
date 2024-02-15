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
        .splitCsv ( header:true, sep:',' )
        .map { create_channel(it) }
        .set { meta }

    emit:
    meta
}

def create_channel(LinkedHashMap row) {
    def meta = [:]
    meta.studyid = "${row.studyid}" == "" ? -1 : "${row.studyid}".toInteger()
    meta.runid = "${row.runid}" == "" ? -1 : "${row.runid}".toInteger()
    // only allow selecting data using the laneid or plexid fields when the runid field
    // is also specified, otherwise it would catch too unspecific datasets.
    // sets default values of -1 for these specific meta fields
    meta.laneid = ((row.runid.toInteger() < 0) || ("${row.laneid}" == "")) ? -1 : "${row.laneid}".toInteger()
    meta.plexid = ((row.runid.toInteger() < 0) || ("${row.plexid}" == "")) ? -1 : "${row.plexid}".toInteger()
    return meta
}