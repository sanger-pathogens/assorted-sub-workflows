process CURATE_CONSENSUS {
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'

    publishDir "${params.outdir}/${meta.id}/curated_consensus", mode: 'copy', overwrite: true

    container 'quay.io/biocontainers/python:3.10.2'

    input:
    tuple val(meta), file(vcf_final), path(reference), path(ref_index)

    output:
    tuple val(meta), path("*.fa"),  emit: curated_consensus

    script:
    ref_basename = reference.baseName
    align_script = "${projectDir}/lib/generate_consensus.py"
    """
    python3 ${align_script} -v '${vcf_final}' -i '${ref_index}' -o '${meta.id}_${ref_basename}.fa' -s '${meta.id}'
    """
}
