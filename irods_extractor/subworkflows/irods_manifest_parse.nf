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
    meta.studyid = row.studyid

    if (row.studyid && row.runid) {
        meta.runid = "${row.runid}"
    }

    if (row.studyid && row.runid && row.laneid) {
        meta.laneid = "${row.laneid}"
    }

    if (row.studyid && row.runid && row.laneid && row.plexid) {
        meta.plexid = "${row.plexid}"
    }

    return meta
}