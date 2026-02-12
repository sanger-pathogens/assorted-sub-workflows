process CHECKM2 {
    tag "${meta.ID}"
    label "cpu_4"
    label "mem_10"
    label "time_queue_from_small"

    container 'quay.io/biocontainers/checkm2:1.0.2--pyh7cba7a3_0'

    input:
    tuple val(meta), path(fastas), path(checkm2_db)

    output:
    tuple val(meta), path(report_txt), emit: results
    tuple val(meta), path(fastas), path(report_txt), emit: results_with_bin

    script:
    report_txt = "${meta.ID}_checkm2_report.tsv"
    """
    checkm2 predict -x .fasta --threads ${task.cpus} --input ${fastas} --database_path ${checkm2_db} --output-directory checkm2

    # move the output file names to something slightly more descriptive
    
    mv checkm2/quality_report.tsv ${report_txt}
    """
}