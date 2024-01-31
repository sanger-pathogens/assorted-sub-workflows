process METADATA {
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'
    
    conda 'anaconda::pandas=2.1.4'
    // NOTE discrepancy of version as v2.1.4 not avialable publicly AFAIK
    container "${ profile.name == 'standard' ? '/software/pathogen/images/python-pandas.simg' : 'biocontainers/pandas:1.5.1_cv1' }"
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