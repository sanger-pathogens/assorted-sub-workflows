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
    gatk --java-options ${params.gatk_java_args} CreateSequenceDictionary \
     --REFERENCE ${reference} \
     --OUTPUT ${reference.baseName}.dict
    """
}

process GATK_HAPLOTYPECALLER {
    label 'cpu_2'
    label 'mem_4'
    label 'time_1'

    container 'broadinstitute/gatk:4.5.0.0'

    publishDir "${params.outdir}/${meta.ID}/gatk", pattern:"*_bamout.*", enabled: params.keep_gatk_bam, mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(sorted_reads), path(sorted_reads_index), path(reference), path(reference_index), path(reference_dict)

    output:
    tuple val(meta), path("${output_vcf}"),  emit: vcf
    tuple val(meta), path("${output_bam}"), path("${output_bai}"),  emit: bamout

    script:
    output_vcf="${meta.ID}.gatk_haplocall_raw.vcf"
    output_bam="${meta.ID}.gatk_haplocall_raw.bam"
    output_bai="${meta.ID}.gatk_haplocall_raw.bai"
    """
    gatk --java-options "${params.gatk_java_args}" HaplotypeCaller  \
      -R ${reference} \
      -I ${sorted_reads} \
      -O ${output_vcf} \
      -bamout ${output_bam} \
      --output-mode ${params.gatk_haplotypecaller_output_mode} \
      --native-pair-hmm-threads ${task.cpus}
    """
}