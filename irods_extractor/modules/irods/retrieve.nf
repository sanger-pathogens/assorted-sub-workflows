process RETRIEVE_CRAM {
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'
    input:
    tuple val(meta), val(cram_path)

    output:
    tuple val(meta), path("*.cram"), emit: path_channel

    script:
    """
    iget -K ${cram_path}
    """
}