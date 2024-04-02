process METADATA {
    label 'cpu_1'
    label 'mem_1'
    label 'time_30m'
    
    conda 'anaconda::pandas=2.1.4'
    // NOTE v2.1.4 not avialable publicly AFAIK so prefering custom image with v2.2.1 vs. quay.io/biocontainers/pandas:1.5.2
    container 'quay.io/sangerpathogens/pandas:2.2.1

    publishDir "${params.outdir}/", mode: 'copy', overwrite: true, pattern: "${timestampout}"

    input:
    path(metadata)

    output:
    path("${timestampout}")
    
    script:
    maptocsv = "${projectDir}/assorted-sub-workflows/irods_extractor/bin/map_to_csv.py"
    timestampout = "${workflow.start}_metadata.csv"
    """
    ${maptocsv} --input_map_list ${metadata}
    mv metadata.csv ${timestampout}
    """
}