process COMBINE_SYLPH_REPORTS {
    tag ""
    label 'cpu_1'
    label 'mem_1'
    label 'time_from_queue_small'

    publishDir "/sylph/", pattern: "sylph_profile.tsv", mode: 'copy', overwrite: true

    container 'ubuntu:22.04'

    input:
    tuple val(meta), path(sylph_reports, stageAs: "reports/*")

    output:
    tuple val(meta), path("sylph_profile.tsv"), emit: sylph_report

    script:
    """
    ls reports/*.tsv | sort | head -n 1 | xargs head -n 1 > sylph_profile.tsv
    tail -q -n +2 reports/*.tsv >> sylph_profile.tsv
    """
}

process GROUP_SYLPH_REFS_BY_TAXON {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_from_queue_small'

    container 'quay.io/sangerpathogens/pandas:2.2.1'

    input:
    tuple val(meta), path(sylph_report), path(sylphtax_report)

    output:
    tuple val(meta), path("${meta.ID}/*.tsv") , emit: taxon_group_ref_reports

    script:
    """
    ${workflow.projectDir}/assorted-sub-workflows/sylph_refset/bin/group_refs.py \\
        --sylph_prof_report ${sylph_report} \\
        --sylphtax_report ${sylphtax_report} \\
        --taxonomic_group ${params.taxonomic_grouping} \\
        --prefix ${meta.ID} \\
        --outdir ${meta.ID}
    """
}

process COMBINE_REFS_ACROSS_SAMPLES {
    tag "${taxon_group}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_from_queue_small'

    publishDir "${params.outdir}/sylph/taxon_group_ref_reports", mode: 'copy', overwrite: true

    container 'ubuntu:22.04'

    input:
    tuple val(taxon_group), path(taxon_group_ref_reports, stageAs: "ref_reports/*")

    output:
    tuple val(taxon_group), path("${taxon_group}.tsv"), emit: taxon_group_ref_report

    script:
    """
    ls ref_reports/*.tsv | head -n 1 | xargs head -n 1 > ${taxon_group}.tsv
    tail -n +2 ref_reports/*.tsv >> ${taxon_group}.tsv
    """
}