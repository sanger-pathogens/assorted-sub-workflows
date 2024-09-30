process DOWNLOAD_METADATA {
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    container 'quay.io/sangerpathogens/enadownloader:v2.3.2-fb2c2cca-bookworm'

    if (params.publish_metadata) {
        publishDir "${params.outdir}/metadata", mode: 'copy', overwrite: true
    }

    input:
    tuple val(meta), path(accessions)

    output:
    tuple val(meta), path("metadata.tsv"), emit: metadata_tsv

    script:
    """
    enadownloader --input ${accessions} --type ${params.accession_type} --write-metadata
    """
}

def write_filter_tsv(filter_map, output_path) {
    // File output = new File(output_path)

    output_path.withWriter('UTF-8') { writer ->
        writer.writeLine("column\tfilter\tdatatype")  // Header

        filter_map.each { col, values ->
            String filter = values[0]
            String datatype = values[1]
            writer.writeLine("${col}\t\"${filter}\"\t${datatype}")
        }
    }
    return output_path
}

process FILTER_METADATA {
    tag "${meta.ID}"
    label "cpu_1"
    label "mem_500M"
    label "time_30m"

    container 'quay.io/sangerpathogens/pysam:0.0.2'

    input:
    tuple val(meta), path(metadata_tsv)
    val(filters_map)
    val(select)
    val(remove_header)

    output:
    tuple val(meta), path(filtered_metadata), emit: filtered_metadata

    script:
    filtered_metadata = "filtered_metadata.tsv"
    filter_tsv = file("filter.tsv")
    write_filter_tsv(filters_map, filter_tsv)
    select_opt = select ? "--select ${select.join(' ')}" : ""
    remove_header_opt = remove_header ? "--remove_header" : ""
    """
    filter_metadata.py \\
        --filter_manifest ${filter_tsv} \\
        --input ${metadata_tsv} \\
        --output ${filtered_metadata} \\
        ${select_opt} \\
        ${remove_header_opt}
    """
}
