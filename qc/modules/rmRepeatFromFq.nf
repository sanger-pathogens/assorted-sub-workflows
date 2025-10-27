// this default script path should work out of the box on any pipeline.
params.script_src_path="${projectDir}/rvi_toolbox/bin/" 

// NOTE: the main reason this param was added is to allow to run tests under
//       rvi_toolbox dir


process RMREPEATFROMFASTQ {
    tag "${meta.id}"
    publishDir enabled: params.debug_preproc_output, "${params.results_dir}/${meta.id}/preprocessing/", mode: "copy", pattern:"*.{trf,fastq}"
    container "quay.io/gsu-pipelines/rvi-vp-basecontainer"

    input:
        tuple val(meta), path(fastq_1), path(fastq_2), path(trf_out_1), path(trf_out_2)

    output:
        tuple val(meta), path("${meta.id}_trf_1.fastq"), path("${meta.id}_trf_2.fastq"), emit: fastqs
        tuple path("combined.trf"), path("unpaired.trf"), emit: combined_trfs

    script:
    """
    ${params.script_src_path}trfCombine.py --trf1 ${trf_out_1} --trf2 ${trf_out_2} -o combined.trf -u unpaired.trf
    ${params.script_src_path}rmRepeatFromFq.py -i ${fastq_1} -t combined.trf -o ${meta.id}_trf_1.fastq
    ${params.script_src_path}rmRepeatFromFq.py -i ${fastq_2} -t combined.trf -o ${meta.id}_trf_2.fastq
    """
}