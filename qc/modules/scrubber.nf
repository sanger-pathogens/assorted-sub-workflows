process SRA_HUMAN_SCRUBBER {
    tag "${meta.ID}"
    label "cpu_1"
    label "mem_10"
    label "time_queue_form_normal"

    container "quay.io/biocontainers/sra-human-scrubber:1.0.2021_05_05--hdfd78af_0"

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