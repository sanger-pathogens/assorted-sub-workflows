def translateKey(in_key) {
    switch(in_key){ 
        case 'studyid':
        return 'study_id'

        case 'runid':
        return 'id_run'

        case 'laneid':
        return 'lane'

        case 'plexid':
        return 'tag_index'

        //if none of the cases return in_key
        default:
        return in_key
    }
}

def avuIdQuery(meta_query) {
    def query_list = []
    // with validation for numeric id types
    meta_query.each { key, value ->
        if (value.class.simpleName == 'String' || value >= 0) {
            def irods_key = translateKey(key)
            def query_part = """{"a": "${irods_key}", "v": "${value}"}"""
            query_list.add(query_part)
         }
    }
    return query_list
}

def objectOrCollection() {
    switch(params.read_type) {
        case 'ont':
            return '"collection": true'
        case 'all':
            return '"object": true, "collection": true'
        case 'illumina':
            return '"object": true'
        default:
            log.error("read_type '${params.read_type}' was not one of illumina|ont|all")
    }
}

process JSON_PREP {
    label 'cpu_1'
    label 'mem_1'
    label 'time_30m'

    container 'quay.io/sangerpathogens/jq:1.6'

    input:
    val(meta)

    output:
    path(json_file), emit: path_channel

    script:
    json_file="input.json"
    """
    jq -n '{op: "metaquery", args: {${objectOrCollection()}, "avu": true}, target: {avus: ${avuIdQuery(meta)}}}' > ${json_file}
    """
}
