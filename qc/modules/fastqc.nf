process FASTQC {
    tag "${meta.ID}"
    label 'cpu_2'
    label 'mem_2'
    label 'time_12'

    publishDir "${params.outdir}/${meta.ID}/fastqqc/", pattern: "*.zip", mode: 'copy', overwrite: true, enabled: params.save_fastqc

    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

    input:
    tuple val(meta), path(read_1), path(read_2)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path(read_1_fastqc), path(read_2_fastqc), emit: zip

    script:
    read_1_fastqc = "${meta.ID}_1_fastqc.zip"
    read_2_fastqc = "${meta.ID}_2_fastqc.zip"
    """
    fastqc \
        -f fastq \
        --threads ${task.cpus} \
        ${read_1} ${read_2}
    """
}
