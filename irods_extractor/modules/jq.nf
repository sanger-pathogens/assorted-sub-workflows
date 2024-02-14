// bellow functions appear to not be thread-safe
// see https://gitlab.internal.sanger.ac.uk/sanger-pathogens/pipelines/assorted-sub-workflows/-/merge_requests/40
// leaving them here if we can find a way to bring them in cleanly
// in the meantime, bring back the less elegant bash function from 7ac9fca2
/*
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
    def query_list = ["""{"a": "target", "v": "1"}, {"a": "type", "v": "cram"}"""]
    // with validation for numeric id types
    meta_query.each { key, value ->
        irods_key = translateKey(key)
        query_part = """{"a": "${irods_key}", "v": "${value}"}"""
        query_list.add(query_part)
        }
    return query_list
}
*/

process JSON_PREP {
    label 'cpu_1'
    label 'mem_1'
    label 'time_30m'

    input:
    val(meta)

    output:
    path(json_file), emit: path_channel

    script:
    json_file="input.json"  
    """
    avu_id_query () {
    avukey="\${1}"
    avuval="\${2}"
    # with validation for numeric id types
    if [ "\${avuval}" -ge 0 ] ; then
            avuq="{a: \\"\${avukey}\\", v: \\"\${avuval}\\"}, "
        else
            avuq=""
        fi
        echo "\${avuq}"
    }
    studyq=\$(avu_id_query 'study_id' ${meta.studyid})
    runq=\$(avu_id_query 'id_run' ${meta.runid})
    laneq=\$(avu_id_query 'lane' ${meta.laneid})
    plexq=\$(avu_id_query 'tag_index' ${meta.plexid})
    jq -n "{op: \\"metaquery\\", args: {object: true, \\"avu\\": true}, target: {avus: [\${studyq}\${runq}\${laneq}\${plexq}{a: \\"target\\", v: \\"1\\"}, {a: \\"type\\", v: \\"cram\\"}]}}" > ${json_file}
    """
}

process JSON_PARSE {
    label 'cpu_1'
    label 'mem_1'
    label 'time_30m'
    
    input:
    path(lane_file)

    output:
    stdout emit: paths
    path("irods_paths.json"), emit: json_file

    script:
    //format is dodgy when it comes off of IRODS so second JQ fixes the formatting
    """
    jq '.result[] | .collection + "/" + .data_object' ${lane_file}
    jq -r '' ${lane_file} > irods_paths.json
    """
}