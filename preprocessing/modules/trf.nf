process TRF {
    tag "${meta.ID}"
    label 'mem_1'
    label 'cpu_1'
    label 'time_1'

    container "quay.io/biocontainers/trf:4.09.1--h031d066_6"

    publishDir enabled: params.publish_intermediate_trf, mode: 'copy', pattern: "*.trf", path: "${params.outdir}/${meta.ID}/preprocessing/intermediates/trf"

    input:
    tuple val(ID), val(read_type), val(meta), path(fastq), path(fasta)

    output:
    tuple val(ID), val(read_type), val(meta), path(fastq), path(trf_marked), emit: fastq_trf_marked

    script:
    trf_marked = "${meta.ID}_${read_type}.trf"
    """
    if [ -f "$fasta" ] && [ ! -s "$fasta" ]; then
        touch "${trf_marked}"
    else
        trf ${fasta} ${params.trf_cli_options} > ${trf_marked}
    fi
    """
}