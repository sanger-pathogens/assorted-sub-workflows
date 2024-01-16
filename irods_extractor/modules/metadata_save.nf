process METADATA {
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'
    
    container '/software/pathogen/images/python-pandas.simg'

    publishDir "${params.results_dir}/", mode: 'copy', overwrite: true, pattern: "metadata.csv"

    input:
    path(metadata)

    output:
    path("metadata.csv")

    script:
    """
    map_to_csv.py --input_map_list ${metadata}
    """
}