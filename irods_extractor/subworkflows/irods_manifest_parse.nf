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
    def studyid = { "${row.studyid}" == "" ? -1 : "${row.studyid}".toInt() }
    meta.studyid = studyid
    def runid = { "${row.runid}" == "" ? -1 : "${row.runid}".toInt() }
    meta.runid = runid
    // only allow selecting data using the laneid or plexid fields when the runid field
    // is also specified, otherwise it would catch too unspecific datasets.
    // sets default values of -1 for these specific meta fields
    def laneid = { ((row.runid < 0) || ("${row.laneid}" == "")) ? -1 : "${row.laneid}".toInt() }
    meta.laneid = laneid
    def plexid = { ((row.runid < 0) || ("${row.plexid}" == "")) ? -1 : "${row.plexid}".toInt() }
    meta.plexid = plexid
    return meta
}