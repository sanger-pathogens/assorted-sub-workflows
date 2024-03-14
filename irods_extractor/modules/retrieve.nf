process RETRIEVE_CRAM {
    label 'cpu_1'
    label 'mem_1'
    label 'time_30m'
    input:
    tuple val(meta), val(cram_path)

    output:
    tuple val(meta), path("*.cram"), emit: path_channel

    // output file renaming logic relies on the iRODS folder structure, which is not perfect, but can't think of anything else there
    script:
    stripcrampath = "${cram_path}".strip().strip('"')
    headcramdir = file(stripcrampath).parent.name
    irodscram = file(stripcrampath).name
    """
    iget -K "${stripcrampath}" && \
    if [[ "${headcramdir}" =~ 'plex' ]] ; then 
      mv ${irodscram} ${headcramdir}_${irodscram}
    fi
    """
}