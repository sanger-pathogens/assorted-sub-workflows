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

def generate_query(params) {

    def keys = ['studyid', 'runid', 'laneid', 'plexid']

    def query_list = keys
        .findAll { params[it] != null }
        .collect { key ->
            def irods_key = translateKey(key)
            "-q ${irods_key}=${params[key]} "
        }

    return query_list.join('')
}


process BEEFEATER {
    label 'cpu_1'
    label 'mem_2'
    label 'time_1'

    container 'quay.io/sangerpathogens/beefeater:v1.0.1-6b5c227e'

    output:
    path("*output.csv"), emit: csv_ch //this is a json but the output file name is messed up to fix

    script:
    def query    = generate_query(params)
    def manifest = params.manifest_of_lanes ? "--manifest ${params.manifest_of_lanes}" : ""
    def search   = params.search ? "" : "--get"
    """
    beefeater search \\
        ${query} \\
        ${manifest} \\
        ${search} \\
        --file_format json
    """
}
