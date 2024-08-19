process RETRIEVE_CRAM {
    label 'cpu_1'
    label 'mem_1'
    label 'time_30m'

    //using the singularity image below the ISG/experiemntal/irods/4.3.0 module
    container '/software/isg/private/experimental/irods/4.3.0/lib/4.3.0_135981_feat-inital.sif'

    input:
    tuple val(meta), val(cram_path)

    output:
    tuple val(meta), path("*.{cram,bam}"), emit: path_channel

    script:
    """
    iget -K ${cram_path}
    """
}