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
    if ((meta.runid < 0) && (("${row.plexid}" != "") || ("${row.laneid}" != ""))) {
        log.warn ("cannot submit an iRODS query where the laneid or plexid parameters are spcified but not the runid, as this query would catch too many file objects.\nThe row '${row.studyid},${row.runid},${row.laneid},${row.plexid}' in the input manifest is ignored")
        return "none" // using same format of empty channel items as in COMBINED_INPUT L39 
    } else {
        meta.laneid = "${row.laneid}" == "" ? -1 : "${row.laneid}".toInteger()
        meta.plexid = "${row.plexid}" == "" ? -1 : "${row.plexid}".toInteger()
        return meta
    }
}