process BATON {
    label 'cpu_1'
    label 'mem_250M'
    label 'time_30m'
    maxForks = 10

    container "ghcr.io/wtsi-npg/ub-16.04-baton-irods-4.2.7:5.0.1"
    
    input:
    path(json_file)

    output:
    path(lane_file), emit: json_file

    script:
    lane_file="info.json"
    """
    baton-do --file ${json_file} --zone seq | jq -r '' > ${lane_file}
    """
}

process BATON_GET {
    label 'cpu_1'
    label 'mem_100M'
    label 'time_30m'
    maxForks = 10

    container "ghcr.io/wtsi-npg/ub-16.04-baton-irods-4.2.7:5.0.1"
    
    input:
    val(meta)

    output:
    tuple val(meta), path("${meta.ID}{,_${meta.subset}}.{cram,bam}"), emit: path_channel

    script:
    """
    jq -n '{collection: "${meta.collection}", data_object: "${meta.data_obj}"}' | baton-get --save
    """
}

process BATON_GET_COLLECTION {
    label 'cpu_1'
    label 'mem_100M'
    label 'time_30m'
    maxForks = 10

    publishDir path: { if ("${params.save_method}" == "nested") "${params.outdir}/${meta.ID}/${params.raw_reads_prefix}pod5/" else "${params.outdir}/pod5/" }, enabled: params.save_fastqs, mode: 'copy', overwrite: true


    container "ghcr.io/wtsi-npg/ub-16.04-baton-irods-4.2.7:5.0.1"
    
    input:
    val(meta)

    output:
    tuple val(meta), path("*.{pod5,fastq.gz}"), emit: path_channel

    script:
    """
    jq -n '{collection: "${meta.irods_path}"}' \\
    | baton-list --contents \\
    | jq -r '.contents[]' \\
    | baton-get --save
    """
}