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
    sylph sketch --threads ${task.cpu} -1 ${read_1} -2 ${read_2} -k ${params.sketch_size} -d paired_sketches
    """
}

process SYLPH_PROFILE {
    tag "${meta.ID}"
    label 'cpu_2'
    label 'mem_2'
    label 'time_from_queue_small'

    publishDir "${params.outdir}/${meta.ID}/sylph/", pattern: "*.tsv", mode: 'copy', overwrite: true, enabled: true

    container 'gitlab-registry.internal.sanger.ac.uk/sanger-pathogens/docker-images/sylph:0.8.1--ha6fb395_0'

    input:
    tuple val(meta), path(sketch)

    output:
    path "${meta.ID}_sylph_profile.tsv", emit: sylph_report

    script:
    """
    sylph profile --threads ${task.cpu} -o ${meta.ID}_sylph_profile.tsv -k ${params.sketch_size} ${sketch} ${params.sylph_db}
    """
}

process SYLPHTAX_TAXPROF {
    tag "${meta.ID}"
    label 'cpu_2'
    label 'mem_2'
    label 'time_from_queue_small'

    publishDir "${params.outdir}/${meta.ID}/sylph/", pattern: "*.sylphmpa", mode: 'copy', overwrite: true, enabled: true

    container 'gitlab-registry.internal.sanger.ac.uk/sanger-pathogens/docker-images/sylph-tax:1.2.0'

    input:
    tuple val(meta), path(sylph_report)

    output:
    path "${meta.ID}_sylphtax_profile.sylphmpa", emit: sylphtax_mpa_report

    script:
    """
    sylph-tax taxprof --threads ${task.cpu} -o ${meta.ID}_sylphtax_profile ${sylph_report} -t ${params.sylphtax_db_tag}
    """
}