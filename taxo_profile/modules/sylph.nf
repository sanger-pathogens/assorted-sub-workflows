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
    tuple val(meta), path("${meta.ID}.paired.sylsp"), emit: sketch

    script:
    """
    sylph sketch -t ${task.cpus} -1 ${read_1} -2 ${read_2} -k ${params.sketch_size} -S ${meta.ID} -d ./
    """
}

process SYLPH_PROFILE {
    tag "${meta.ID}"
    label 'cpu_2'
    label 'mem_20'
    label 'time_from_queue_small'

    publishDir "${params.outdir}/${meta.ID}/sylph/", pattern: "*.tsv", mode: 'copy', overwrite: true

    container 'gitlab-registry.internal.sanger.ac.uk/sanger-pathogens/docker-images/sylph:0.8.1--ha6fb395_0'

    input:
    tuple val(meta), path(sketch)

    output:
    tuple val(meta), path("${meta.ID}_sylph_profile.tsv"), emit: sylph_report

    script:
    def estimate_unknown = params.sylph_estimate_unknown ? (params.sylph_read_seq_id ? "-u --read-seq-id ${params.sylph_read_seq_id}" : "-u") : ""
    """
    sylph profile -t ${task.cpus} -o ${meta.ID}_sylph_profile.tsv -k ${params.sketch_size} ${sketch} ${params.sylph_db} ${estimate_unknown}
    """
}

process SYLPH_QUERY {
    tag "${meta.ID}"
    label 'cpu_2'
    label 'mem_20'
    label 'time_from_queue_small'

    publishDir "${params.outdir}/${meta.ID}/sylph/", pattern: "*.tsv", mode: 'copy', overwrite: true

    container 'gitlab-registry.internal.sanger.ac.uk/sanger-pathogens/docker-images/sylph:0.8.1--ha6fb395_0'

    input:
    tuple val(meta), path(sketch)

    output:
    tuple val(meta), path("${meta.ID}_sylph_profile.tsv"), emit: sylph_report

    script:
    """
    sylph query ${sketch} ${params.sylph_db} -t ${task.cpus} -o ${meta.ID}_sylph_profile.tsv
    """
}

process SYLPHTAX_TAXPROF {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_2'
    label 'time_from_queue_small'

    publishDir "${params.outdir}/${meta.ID}/sylph/", pattern: "*.sylphmpa", mode: 'copy', overwrite: true

    container 'gitlab-registry.internal.sanger.ac.uk/sanger-pathogens/docker-images/sylph-tax:1.2.0'

    input:
    tuple val(meta), path(sylph_report)

    output:
    tuple val(meta), path("${meta.ID}_sylphtax_profile.sylphmpa") , emit: sylphtax_mpa_report

    script:
    """
    sylph-tax taxprof ${sylph_report} -t ${params.sylphtax_db_tag}
    mv ${meta.ID}.sylphmpa ${meta.ID}_sylphtax_profile.sylphmpa
    """
}

process SYLPH_SUMMARIZE {
    label 'cpu_1'
    label 'mem_4'
    label 'time_queue_from_small'

    container 'quay.io/sangerpathogens/pandas:2.2.1'

    publishDir "${params.outdir}/sylph/", mode: 'copy', overwrite: true

    input:
    path(sylph_reports)

    output:
    path("references.txt"), emit: references
    path("sylph_summary.tsv"), emit: sylph_summary

    script:
    // Filter once with thresholds.
    """
    ${workflow.projectDir}/assorted-sub-workflows/sylph_db_refinement/bin/sylph_summarize.py \
        --reports ${sylph_reports} \
        --ani ${params.sylph_ani} \
        --cov ${params.sylph_cov} \
        --ani-column Adjusted_ANI \
        --cov-column Eff_cov \
        --out-references references.txt \
        --out-summary sylph_summary.tsv
    """
}
