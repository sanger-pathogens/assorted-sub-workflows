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

    // see irods.config params scope for recommended default values to set at pipeline level
    switch(meta.size()) {
        case 1:
            """
            jq -n '{op: "metaquery", args: {object: true, "avu": true}, target: {avus: [{a: "study_id", v: "${meta.study}"}, {a: "target", v: "1"}, {a: "type", v: "cram"}]}}' > ${json_file}
            """
            break
        case 2:
            """
            jq -n '{op: "metaquery", args: {object: true, "avu": true}, target: {avus: [{a: "study_id", v: "${meta.study}"}, {a: "id_run", v: "${meta.runid}"}, {a: "target", v: "1"}, {a: "type", v: "cram"}]}}' > ${json_file}
            """
            break
        case 3:
            """
            jq -n '{op: "metaquery", args: {object: true, "avu": true}, target: {avus: [{a: "study_id", v: "${meta.study}"}, {a: "id_run", v: "${meta.runid}"}, {a: "lane", v: "${meta.laneid}"}, {a: "target", v: "1"}, {a: "type", v: "cram"}]}}' > ${json_file}
            """
            break
        case 4:
            """
            jq -n '{op: "metaquery", args: {object: true, "avu": true}, target: {avus: [{a: "study_id", v: "${meta.study}"}, {a: "id_run", v: "${meta.runid}"}, {a: "lane", v: "${meta.laneid}"}, {a: "tag_index", v: "${meta.plexid}"}, {a: "target", v: "1"}, {a: "type", v: "cram"}]}}' > ${json_file}
            """
            break
    }
}

process JSON_PARSE {
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    publishDir "${params.outdir}/", mode: 'copy', overwrite: true, pattern: "irods_paths.json"

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