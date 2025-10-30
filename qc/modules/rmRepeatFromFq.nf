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
    tuple val(meta), path(output_1), path(output_2), emit: fastqs
    tuple val(meta), path(output_1_gz), path(output_2_gz), emit: fastqs_gz   
    tuple path("combined.trf"), path("unpaired.trf"), emit: combined_trfs

    script:
    output_1="${meta.ID}_trf_1.fastq"
    output_2="${meta.ID}_trf_2.fastq"
    output_1_gz="${output_1}.gz"
    output_2_gz="${output_2}.gz"
    """
    ${params.script_src_path}trfCombine.py --trf1 ${trf_out_1} --trf2 ${trf_out_2} -o combined.trf -u unpaired.trf
    ${params.script_src_path}rmRepeatFromFq.py -i ${fastq_1} -t combined.trf -o ${output_1}
    ${params.script_src_path}rmRepeatFromFq.py -i ${fastq_2} -t combined.trf -o ${output_2}
    gzip -c "${output_1}" > ${output_1}.tmp.gz
    gzip -c "${output_2}" > ${output_2}.tmp.gz
    mv ${output_1}.tmp.gz ${output_1_gz}
    mv ${output_2}.tmp.gz ${output_2_gz}
    """
}