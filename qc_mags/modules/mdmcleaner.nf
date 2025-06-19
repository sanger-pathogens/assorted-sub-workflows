process MDMCLEANER {
    tag "${meta.ID}"
    label "cpu_4"
    label "mem_8"
    label "time_12"

    container  'quay.io/biocontainers/mdmcleaner:0.8.7--pyh7cba7a3_0'

    publishDir mode: 'copy', path: "${params.outdir}/gunc/"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path(mdmcleaner_output), emit: results

    script:
    """
    mdmcleaner set_configs --db_basedir ${params.mdmcleaner_db} --threads ${task.cpus}
    mdmcleaner clean -i ${fasta}/* -o mdmcleaner_output -t ${task.cpus} -c mdmcleaner.config
    """
}