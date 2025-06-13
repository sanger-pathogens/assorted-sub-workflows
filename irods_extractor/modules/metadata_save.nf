process METADATA {
    label 'cpu_1'
    label 'mem_1'
    label 'time_30m'
    
    conda 'anaconda::pandas=2.2.1'
    container 'quay.io/sangerpathogens/pandas:2.2.1'

    publishDir "${params.outdir}/", mode: 'copy', overwrite: true

    input:
    path(metadata)
    val(metadata_tag)

    output:
    path("metadata_${metadata_tag}_*.csv")
    
    script:
    def maptocsv = "${projectDir}/assorted-sub-workflows/irods_extractor/bin/map_to_csv.py"
    def dateFormat = params.short_metacsv_name ? "+%Y-%m-%d" : "+%Y-%m-%dT%H-%M-%S"
    """
    DATE=\$(date ${dateFormat})
    OUTPUT_FILE="metadata_${metadata_tag}_\${DATE}.csv"

    ${maptocsv} --input_map_list ${metadata} --output "\$OUTPUT_FILE"
    """
}
