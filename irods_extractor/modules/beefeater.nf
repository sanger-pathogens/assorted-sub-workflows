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

def query(meta_query) {
    query_list = []
    meta_query.each { key, value ->
        if (value.class.simpleName == 'String' || value >= 0){
            def irods_key = translateKey(key)
            def query_part = "-q ${irods_key}=${value}"
            query_list.add(query_part)
         }
    }
    return query_list.join(' ') 
}


process BEEFEATER {
    label 'cpu_1'
    label 'mem_10'
    label 'time_30m'

    container 'quay.io/sangerpathogens/beefeater:0.0.1-c1'

    input:
    val(meta)

    output:
    path("*output.json"), emit: csv_ch

    script:
    def search_option = params.search ? "" : "--get"

    """
    beefeater search ${query(meta)} ${search_option} --file_format json
    """
}
