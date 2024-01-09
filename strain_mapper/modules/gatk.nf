process GATK_REF_DICT {
    label 'cpu_2'
    label 'mem_4'
    label 'time_1'

    container 'broadinstitute/gatk:4.5.0.0'

    publishDir "${params.outdir}/sorted_ref", mode: 'copy', overwrite: true, pattern: '*.dict'

    input:
    path(reference)

    output:
    path("*.dict"), emit: ref_dict

    script:
    """
    gatk --java-options "${params.gatk_java_args}" CreateSequenceDictionary \
     --REFERENCE ${reference} \
     --OUTPUT ${reference.baseName}.dict
    """
}

process GATK_HAPLOTYPECALLER {
    label 'cpu_2'
    label 'mem_4'
    label 'time_1'

    container 'broadinstitute/gatk:4.5.0.0'

    publishDir "${params.outdir}/${meta.id}/gatk", pattern:"*_bamout.bam", enabled: params.keep_gatk_bam, mode: 'copy', overwrite: true
    publishDir "${params.outdir}/${meta.id}/gatk", pattern:"*.vcf.gz", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(sorted_reads), path(sorted_reads_index), path(reference), path(reference_index), path(reference_dict)

    output:
    tuple val(meta), path("${output_vcf}"),  emit: vcf
    tuple val(meta), path("${bamout}"),  emit: bamout

    script:
    output_vcf="${meta.id}.vcf.gz"
    bamout="${meta.id}_bamout.bam"
    """
    gatk --java-options "${params.gatk_java_args}" HaplotypeCaller  \
      -R ${reference} \
      -I ${sorted_reads} \
      -O ${output_vcf} \
      -bamout ${bamout} \
      --native-pair-hmm-threads ${task.cpus}
    """
}