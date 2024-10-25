process METADATA {
    label 'cpu_1'
    label 'mem_1'
    label 'time_30m'
    
    conda 'anaconda::pandas=2.2.1'
    container 'quay.io/sangerpathogens/pandas:2.2.1'

    publishDir "${params.outdir}/", mode: 'copy', overwrite: true, pattern: "${timestampout}"

    input:
    path(metadata)
    val(metadata_tag)

    output:
    path("${timestampout}")
    
    script:
    maptocsv = "${projectDir}/assorted-sub-workflows/irods_extractor/bin/map_to_csv.py"
    date = "${workflow.start}".split('.')[0] //I would love to use short date not time split on T but will overwrite too much . is miliseconds
    timestampout = "metadata_${metadata_tag}_${date}.csv"
    """
    ${maptocsv} --input_map_list ${metadata} --output "${timestampout}"
    """
}
