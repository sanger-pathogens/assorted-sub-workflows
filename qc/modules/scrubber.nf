process SRA_HUMAN_SCRUBBER {
    tag "${meta.id}"
    label "sra_human_scrubber"

    container "quay.io/gsu-pipelines/rvi-pp-sra-human-scrubber:v1.0"
    input:
        tuple val(meta), path(interleaved_fastq)

    output:
        tuple val(meta), path("${meta.id}_interleaved_clean.fastq") //path("${meta.id}_1_clean.fastq"), path("${meta.id}_2_clean.fastq")

    script:
    """
    scrub.sh -x -s -p ${task.cpus} \
      -i ${interleaved_fastq} -o ${meta.id}_interleaved_clean.fastq
    """
}