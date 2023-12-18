//
// Check input manifest and assign idenifier channels
// Expected to be used as submodule to input data to IRODS extractor
//

workflow IRODS_MANIFEST_PARSE {

    main:
    Channel
        .fromPath( params.manifest )
        .ifEmpty {exit 1, "Cannot find path file ${params.manifest}"}
        .splitCsv ( header:true, sep:',' )
        .map { create_channel(it) }
        .set { meta }

    emit:
    meta
}

def create_channel(LinkedHashMap row) {
    def meta = [:]
    meta.study = row.study

    if (row.study && row.runid) {
        meta.runid = "${row.runid}"
    }

    if (row.study && row.runid && row.laneid) {
        meta.laneid = "${row.laneid}"
    }

    if (row.study && row.runid && row.laneid && row.plexid) {
        meta.plexid = "${row.plexid}"
    }

    return meta
}