process JSON_PREP {
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    input:
    val(study)

    output:
    path(json_file), emit: path_channel

    script:
    json_file="input.json"
    if (study.isNumber())
        """
        jq -n '{op: "metaquery", args: {object: true}, target: {avus: [{a: "study_id", v: "${study}"}, {a: "target", v: "1"}, {a: "type", v: "cram"}]}}' > ${json_file}
        """
    else if (!study.isNumber())
        """
        jq -n '{op: "metaquery", args: {object: true}, target: {avus: [{a: "study", v: "${study}"}, {a: "target", v: "1"}, {a: "type", v: "cram"}]}}' > ${json_file}
        """
    else
        error "unrecognised study input"
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