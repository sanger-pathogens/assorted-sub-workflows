process METADATA {
    label 'cpu_1'
    label 'mem_1'
    label 'time_30m'
    
    conda 'anaconda::pandas=2.1.4'
    // NOTE v2.1.4 not avialable publicly AFAIK so prefering custom image on farm with v2.1.4 vs. quay.io/biocontainers/pandas:1.5.2
    container '/software/pathogen/images/python-pandas.simg'

    publishDir "${params.outdir}/", mode: 'copy', overwrite: true, pattern: "metadata.csv"

    input:
    path(metadata)

    output:
    path("metadata.csv")
    
    script:
    maptocsv = "${projectDir}/assorted-sub-workflows/irods_extractor/bin/map_to_csv.py"
    """
    ${maptocsv} --input_map_list ${metadata}
    """
}