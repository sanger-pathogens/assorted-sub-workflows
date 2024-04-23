process CURATE_CONSENSUS {
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'

    publishDir "${params.outdir}/${meta.ID}/curated_consensus", mode: 'copy', overwrite: true

    conda 'conda-forge::python=3.10.2'
    container 'quay.io/sangerpathogens/pysam:0.0.2'

    input:
    tuple val(meta), file(vcf_final), path(reference), path(ref_index)

    output:
    tuple val(meta), path("*.fa"),  emit: curated_consensus

    script:
    ref_basename = reference.baseName
    align_script = "${projectDir}/assorted-sub-workflows/strain_mapper/bin/vcf_to_psuedo.py"
    """
    python3 ${align_script} -b '${vcf_final}' -r '${reference}' -o '${meta.ID}_${ref_basename}.fa'
    """
}
