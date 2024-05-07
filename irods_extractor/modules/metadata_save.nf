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
    timestampout = "metadata_${metadata_tag}_${workflow.start}.csv"
    """
    ${maptocsv} --input_map_list ${metadata}
    mv metadata.csv ${timestampout}
    """
}
