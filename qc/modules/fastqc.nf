process FASTQC {
    tag "${meta.ID}"
    label 'cpu_2'
    label 'mem_2'
    label 'time_12'

    publishDir "${params.outdir}/${meta.ID}/fastqc/", pattern: "*.zip", mode: 'copy', overwrite: true, enabled: params.save_fastqc

    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

    input:
    tuple val(meta), path(read_1), path(read_2)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*_1_fastqc.zip"), path("*_2_fastqc.zip"), emit: zip

    script:
    """
    fastqc \
        -f fastq \
        --threads ${task.cpus} \
        ${read_1} ${read_2}
    """
}
