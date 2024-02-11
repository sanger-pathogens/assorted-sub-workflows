process JSON_PREP {
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    input:
    val(meta)

    output:
    path(json_file), emit: path_channel

    script:
    json_file="input.json"

    def avu_id_query(avukey, avuval){
        // with validation for numeric id types
        if (avuval > 0) {
            avuq = """{a: "${avukey}", v: "${avuval}"}, """
        else
            avuq = ""
        }
        return avuq
    }
    studyq = avu_id_query('study_id', meta.studyid)
    runq = avu_id_query('id_run', meta.runid)
    laneq = avu_id_query('lane', meta.laneid)
    plexq = avu_id_query('tag_index', meta.plexid)
    // would be brilliant if could iterate over keys in meta without having to know them in advance

    """
    jq -n '{op: "metaquery", args: {object: true, "avu": true}, target: {avus: [${studyq}${runq}${laneq}${plexq}, {a: "target", v: "1"}, {a: "type", v: "cram"}]}}' > ${json_file}
    """
}

process JSON_PARSE {
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'
    
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