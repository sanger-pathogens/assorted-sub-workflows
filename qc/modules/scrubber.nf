process SRA_HUMAN_SCRUBBER {
    tag "${meta.ID}"
    label "cpu_1"
    label "mem_10"
    label "time_queue_form_normal"


    container "quay.io/gsu-pipelines/rvi-pp-sra-human-scrubber:v1.0"
    input:
        tuple val(meta), path(interleaved_fastq)

    output:
        tuple val(meta), path("${meta.ID}_interleaved_clean.fastq") //path("${meta.ID}_1_clean.fastq"), path("${meta.ID}_2_clean.fastq")

    script:
    """
    scrub.sh -x -s -p ${task.cpus} \
      -i ${interleaved_fastq} -o ${meta.ID}_interleaved_clean.fastq
    """
}