process METADATA {
    label 'cpu_1'
    label 'mem_1'
    label 'time_30m'
    
    conda 'anaconda::pandas=2.2.1'
    container 'quay.io/sangerpathogens/pandas:2.2.1'

    publishDir "${params.outdir}/", mode: 'copy', overwrite: false, pattern: "${timestampout}"

    input:
    path(metadata)
    val(metadata_tag)

    output:
    path(timestampout)
    
    script:
    def maptocsv = "${projectDir}/assorted-sub-workflows/irods_extractor/bin/map_to_csv.py"
    def date = params.short_metacsv_name ? "${workflow.start}".split('T')[0] : "${workflow.start}"
    timestampout = "metadata_${metadata_tag}_${date}.csv"
    """
    ${maptocsv} --input_map_list ${metadata} --output "${timestampout}"
    """
}
