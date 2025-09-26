process SYLPH_SKETCH {
    tag "${meta.ID}"
    label 'cpu_2'
    label 'mem_2'
    label 'time_from_queue_small'

    publishDir "${params.outdir}/${meta.ID}/sylph/", pattern: "*.sylsp", mode: 'copy', overwrite: true, enabled: params.save_sylph_sketches

    container 'gitlab-registry.internal.sanger.ac.uk/sanger-pathogens/docker-images/sylph:0.8.1--ha6fb395_0'

    input:
    tuple val(meta), path(read_1), path(read_2)

    output:
    tuple val(meta), path("paired_sketches/${meta.ID}.sylsp"), emit: sketch

    script:
    reads1 = ${meta.ID}.1.fq
    reads2 = ${meta.ID}.2.fq

    """
    sylph sketch --threads ${task.cpu} -1 ${read_1} -2 ${read_2} -d paired_sketches
    """
}

process SYLPH_PROFILE {
    tag "${meta.ID}"
    label 'cpu_2'
    label 'mem_2'
    label 'time_from_queue_small'

    publishDir "${params.outdir}/${meta.ID}/sylph/", pattern: "*.zip", mode: 'copy', overwrite: true, enabled: true

    container 'gitlab-registry.internal.sanger.ac.uk/sanger-pathogens/docker-images/sylph:0.8.1--ha6fb395_0'

    input:
    tuple val(meta), path(sketch)

    output:
    path "${meta.ID}_sylph_profile.tsv", emit: sylph_report

    script:
    """
    sylph profile --threads ${task.cpu} -o ${meta.ID}_sylph_profile.tsv ${sketch} ${params.sylph_db}
    """
}