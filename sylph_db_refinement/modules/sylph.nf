process SYLPH_SKETCH {
    label 'cpu_2'
    label 'mem_16'
    label 'time_from_queue_normal'

    publishDir "${params.outdir}/sylph/", pattern: "*.sylsp", mode: 'copy', overwrite: true, enabled: params.save_sylph_sketches

    container 'gitlab-registry.internal.sanger.ac.uk/sanger-pathogens/docker-images/sylph:0.8.1--ha6fb395_0'

    input:
    tuple val(meta), path(read1), path(read2)

    output:
    tuple val(meta), path("${meta.ID}.sylsp"), emit: sylsp

    script:
    """
    sylph sketch -t ${task.cpus} -1 ${read1} -2 ${read2} -S ${meta.ID} -d ./ -k ${params.sketch_size}
    """
}

// Sylph Version >=0.6 you can do direct profiling of paired-end reads without sketching.
// query: uses ANI and k-mers to determine which species are in the sample.

process SYLPH_QUERY {
    label 'cpu_2'
    label 'mem_20'
    label 'time_from_queue_small'

    publishDir "${params.outdir}/sylph/", pattern: "*_sylph_profile.tsv", mode: 'copy', overwrite: true

    container 'gitlab-registry.internal.sanger.ac.uk/sanger-pathogens/docker-images/sylph:0.8.1--ha6fb395_0'

    input:
    tuple val(meta), path(sample_sketch)
    path(sylph_db)

    output:
    tuple val(meta), path("${meta.ID}_sylph_profile.tsv"), emit: sylph_report

    script:
    // Query a sample sketch against a single chosen .syldb (custom or default).
    """
    sylph query ${sample_sketch} ${sylph_db} -t ${task.cpus} -o ${meta.ID}_sylph_profile.tsv
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
