process CHECKM {
    tag "${meta.ID}"
    label "cpu_4"
    label "mem_64"
    label "time_1"

    container 'quay.io/biocontainers/checkm-genome:1.2.4--pyhdfd78af_0'

    input:
    tuple val(meta), path(fastas)

    output:
    tuple val(meta), path(fastas), path(report_txt), emit: results

    script:
    report_txt = "${meta.ID}/storage/bin_stats_ext.tsv"
    """
    mkdir tmp
    checkm lineage_wf -x fasta ${fastas} ${meta.ID} -t ${task.cpus} --tmpdir tmp --pplacer_threads ${task.cpus}
    """
}

process SUMMARISE_CHECKM {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_250M'
    label 'time_30m'

    container 'quay.io/sangerpathogens/python-curl:3.11'

    input:
    tuple val(meta), path(fasta), path(report)

    output:
    tuple val(meta), path(fasta), path(summary), emit: merged_bins

    script:
    command = "${projectDir}/assorted-sub-workflows/generate_mags/bin_refinement/bin/summarise_checkm.py"
    summary = "${meta.ID}_checkm_summary.tsv"
    """
    ${command} ${report} _ > ${summary}
    """    
}