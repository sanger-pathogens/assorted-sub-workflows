process CHECKM2 {
    tag "${meta.ID}"
    label "cpu_4"
    label "mem_8"
    label "time_30m"

    container 'quay.io/biocontainers/checkm2:1.0.2--pyh7cba7a3_0'

    input:
    tuple val(meta), path(fastas)

    output:
    tuple val(meta), path(report_txt), emit: results
    tuple val(meta), path(fastas), path(report_txt), emit: results_with_bin

    script:
    report_txt = "${meta.ID}_${bin_name}_checkm2_report.tsv"
    """
    checkm2 predict -x .fasta --threads ${task.cpus} --input ${fastas} --output-directory checkm2 --database_path ${params.checkm2_db}

    # move the output file names to something slightly more descriptive
    
    mv checkm2/quality_report.tsv ${report_txt}
    """
}