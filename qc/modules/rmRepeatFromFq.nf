process RMREPEATFROMFASTQ {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    publishDir enabled: params.debug_preproc_output, "${params.outdir}/${meta.ID}/preprocessing/", mode: "copy", pattern:"*.{trf,fastq}"
    container "quay.io/gsu-pipelines/rvi-vp-basecontainer"

    input:
    tuple val(meta), path(fastq_1), path(fastq_2), path(trf_out_1), path(trf_out_2)

    output:
    tuple val(meta), path("${meta.ID}_trf_1.fastq"), path("${meta.ID}_trf_2.fastq"), emit: fastqs
    tuple path("combined.trf"), path("unpaired.trf"), emit: combined_trfs

    script:
    """
    ${params.script_src_path}trfCombine.py --trf1 ${trf_out_1} --trf2 ${trf_out_2} -o combined.trf -u unpaired.trf
    ${params.script_src_path}rmRepeatFromFq.py -i ${fastq_1} -t combined.trf -o ${meta.ID}_trf_1.fastq
    ${params.script_src_path}rmRepeatFromFq.py -i ${fastq_2} -t combined.trf -o ${meta.ID}_trf_2.fastq
    """
}