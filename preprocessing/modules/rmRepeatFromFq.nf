process RM_REPEAT_FROM_PAIRED {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    publishDir enabled: params.publish_intermediate_trf, "${params.outdir}/${meta.ID}/preprocessing/", mode: "copy", pattern:"*.{trf,fastq}"
    container 'quay.io/sangerpathogens/python-curl:3.11'

    input:
    tuple val(meta), path(fastq_1), path(fastq_2), path(trf_out_1), path(trf_out_2)

    output:
    tuple val(meta), path(output_1), path(output_2), emit: fastqs
    tuple path("combined.trf"), path("unpaired.trf"), emit: combined_trfs

    script:
    output_1="${meta.ID}_trf_1.fastq"
    output_2="${meta.ID}_trf_2.fastq"
    """
    ${params.script_src_path}trfCombine.py --trf1 ${trf_out_1} --trf2 ${trf_out_2} -o combined.trf -u unpaired.trf
    ${params.script_src_path}rmRepeatFromFq.py -i ${fastq_1} -t combined.trf -o ${output_1}
    ${params.script_src_path}rmRepeatFromFq.py -i ${fastq_2} -t combined.trf -o ${output_2}
    """
}

process RM_REPEAT_FROM_UNPAIRED {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    publishDir enabled: params.publish_intermediate_trf, "${params.outdir}/${meta.ID}/preprocessing/", mode: "copy", pattern:"*.{trf,fastq}"
    container 'quay.io/sangerpathogens/python-curl:3.11'

    input:
    tuple val(meta), path(fastq), path(trf_out)

    output:
    tuple val(meta), path(output_unpaired), emit: fastqs
    tuple val(meta), path("unpaired_unpaired.trf"), emit: trfs

    script:
    output_unpaired="${meta.ID}_trf_unpaired.fastq"
    """
    ${params.script_src_path}rmRepeatFromFq.py -i ${fastq} -t ${trf_out} -o ${output_unpaired}

    mv ${trf_out} unpaired_unpaired.trf
    """
}