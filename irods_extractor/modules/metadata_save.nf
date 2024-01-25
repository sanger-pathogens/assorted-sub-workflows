process METADATA {
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'
    
    container '/software/pathogen/images/python-pandas.simg'

    publishDir "${params.outdir}/", mode: 'copy', overwrite: true, pattern: "metadata.csv"

    input:
    path(metadata)

    output:
    path("metadata.csv")
    
    maptocsv = "${projectDir}/assorted-sub-workflows/irods_extractor/bin/map_to_csv.py"
    script:
    """
    ${maptocsv} --input_map_list ${metadata}
    """
}