process CURATE_CONSENSUS {
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'

    publishDir "${params.outdir}/${meta.ID}/curated_consensus", mode: 'copy', overwrite: true

    conda 'conda-forge::python=3.10.2'
    container 'quay.io/biocontainers/python:3.10.2'

    input:
    tuple val(meta), file(vcf_final), path(reference), path(ref_index)

    output:
    tuple val(meta), path("*.fa"),  emit: curated_consensus
    tuple val(meta), val("workflow_finished"), emit: finished_ch

    script:
    ref_basename = reference.baseName
    align_script = "${projectDir}/assorted-sub-workflows/strain_mapper/bin/generate_consensus.py"
    """
    python3 ${align_script} -v '${vcf_final}' -i '${ref_index}' -o '${meta.ID}_${ref_basename}.fa' -s '${meta.ID}'
    """
}
