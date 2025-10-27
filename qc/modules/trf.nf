// those are the default option used by kneaddata
params.trf_cli_options= '2 7 7 80 10 50 500 -h -ngs' //"2 5 7 80 10 50 2000 -h -ngs"

process TRF {
    tag "${meta.id}"
    label 'mem_1'
    label 'time_1'

    container "quay.io/biocontainers/trf:4.09.1--h031d066_6"

    publishDir enabled: params.debug_preproc_output, mode: 'copy', failOnError: true, pattern: "*.trf", path: "${params.results_dir}/${meta.id}/preprocessing/"

    input:
    tuple val(meta), path(fasta_R1), path(fasta_R2)

    output:
    tuple val(meta), path("${meta.id}_1.trf"), path("${meta.id}_2.trf"), emit: paired_trf

    script:
    """
    # Handle R1 file
    if [ -f "$fasta_R1" ] && [ ! -s "$fasta_R1" ]; then
        echo "R1 file exists but is empty, writing empty output..."
        touch "${meta.id}_1.trf"
    else
        trf ${fasta_R1} ${params.trf_cli_options} > ${meta.id}_1.trf
    fi

    # Handle R2 file
    if [ -f "$fasta_R2" ] && [ ! -s "$fasta_R2" ]; then
        echo "R2 file exists but is empty, writing empty output..."
        touch "${meta.id}_2.trf"
    else
        trf ${fasta_R2} ${params.trf_cli_options} > ${meta.id}_2.trf
    fi
    """
}