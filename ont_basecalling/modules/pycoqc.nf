process PYCOQC {
    tag "${meta.ID}"

    label 'cpu_1'
    label 'mem_500M'
    label 'time_30m'

    container "quay.io/biocontainers/pycoqc:2.5.2--py_0"
    publishDir "${params.outdir}/qc/${meta.ID}/", mode: 'copy', overwrite: true
    
    input:
    tuple val(meta), path(sequence_summary_file)

    output:
    path("*.html"), emit: html
    path("*.json"), emit: json
        
    script:
    final_qc_file = "${meta.ID}_summary_pycoqc"
    """
    pycoQC -f ${sequence_summary_file} -o ${final_qc_file}.html -j ${final_qc_file}.json
    """
}