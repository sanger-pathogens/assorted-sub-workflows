process FILTER_METADATA {
    tag "${meta.ID}"
    label "cpu_1"
    label "mem_500M"
    label "time_30m"

    container 'quay.io/sangerpathogens/pysam:0.0.2'

    input:
    tuple val(meta), path(metadata_tsv)
    path(filter_manifest)
    val(select)
    val(remove_header)
    val(drop_duplicates_from)

    output:
    tuple val(meta), path(filtered_metadata), emit: filtered_metadata

    script:
    filtered_metadata = "filtered_metadata.tsv"
    select_opt = select ? "--select ${select.join(' ')}" : ""
    remove_header_opt = remove_header ? "--remove_header" : ""
    drop_duplicates_opt = drop_duplicates_from ? "--drop_duplicates_from ${drop_duplicates_from.join(' ')}" : ""

    filter_script = "${projectDir}/assorted-sub-workflows/combined_input/bin/filter_metadata.py"

    """
    ${filter_script} \\
        --filter_manifest ${filter_manifest} \\
        --input ${metadata_tsv} \\
        --output ${filtered_metadata} \\
        ${select_opt} \\
        ${remove_header_opt} \\
        ${drop_duplicates_opt}
    """
}
