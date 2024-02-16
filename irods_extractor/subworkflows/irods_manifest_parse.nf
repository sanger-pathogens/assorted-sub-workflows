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

    //last entry on the row manifest if empty comes in as null so just remove it as its annoying
    row.each { key, value ->
        if (value == null) {
            row[key] = ""
        }
    }

    /*
    these are easy just replace with -1 if we don't want em
    although looking at it why not just drop them if not wanted however not going to change that yet :)
    */

    meta.studyid = "${row.studyid}" == "" ? -1 : "${row.studyid}".toInteger()
    meta.runid = "${row.runid}" == "" ? -1 : "${row.runid}".toInteger()

    /* 
    only allow selecting data using the laneid or plexid fields when the runid field
    is also specified, otherwise it would catch too unspecific datasets.
    sets default values of -1 for these specific meta fields
    */

    def runidInteger = row.runid.isNumber() ? row.runid.toInteger() : -1

    meta.laneid = (runidInteger > 0 && "${row.laneid}".trim() == "") ? -1 : "${row.laneid}".toInteger()

    meta.plexid = (runidInteger > 0 && "${row.plexid}".trim() == "") ? -1 : "${row.plexid}".toInteger()

    return meta
}