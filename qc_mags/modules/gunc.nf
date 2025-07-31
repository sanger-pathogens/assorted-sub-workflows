process GUNC {
    tag "${meta.ID}"
    label "cpu_8"
    label "mem_8"
    label "time_12"

    container  'quay.io/biocontainers/gunc:1.0.6--pyhdfd78af_0'

    publishDir mode: 'copy', path: "${params.outdir}/gunc/"

    input:
    tuple val(meta), path(fastas, stageAs: "fastas/*")

    output:
    tuple val(meta), path(report_tsv), emit: results

    script:
    report_tsv = "${meta.ID}_gunc.tsv"
    """
    mkdir ${meta.ID}_gunc
    gunc run -r ${params.gunc_db} -d fastas -o ${meta.ID}_gunc -t ${task.cpus}

    mv ${meta.ID}_gunc/GUNC.gtdb_95.maxCSS_level.tsv ${report_tsv}
    """
}