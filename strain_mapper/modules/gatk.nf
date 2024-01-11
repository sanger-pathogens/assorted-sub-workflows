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

    publishDir "${params.outdir}/${meta.id}/gatk", pattern:"*_bamout.*", enabled: params.keep_gatk_bam, mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(sorted_reads), path(sorted_reads_index), path(reference), path(reference_index), path(reference_dict)

    output:
    tuple val(meta), path("${output_vcf}"),  emit: vcf
    tuple val(meta), path("${output_bam}"), path("${output_bai}"),  emit: bamout

    script:
    output_vcf="${meta.id}.vcf"
    output_bam="${meta.id}_bamout.bam"
    output_bai="${meta.id}_bamout.bai"
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

process GATK_FILTERING {
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'

    container 'quay.io/biocontainers/bcftools:1.16--haef29d1_2'

    input:
    tuple val(meta), file(vcf_allpos)

    output:
    tuple val(meta), path("${filtered_vcf_allpos}"),  emit: filtered_vcf_allpos

    script:
    filtered_vcf_allpos = "${meta.id}_filtered.vcf"
    """
    bcftools view -o ${filtered_vcf_allpos} \
                  -O 'v' \
                  -i '${params.gatk_vcf_filter}' \
                  '${vcf_allpos}'
    """
}

// TODO The below functions are effectively duplicate code, but with a different publishDir, similar to BCFTOOLS_RAW_VCF, BCFTOOLS_FINAL_VCF...
// That's a lot of duplication...
process GATK_RAW_VCF {
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'
    
    publishDir "${params.outdir}/${meta.id}/gatk/raw_vcf", mode: 'copy', overwrite: true

    container 'quay.io/biocontainers/bcftools:1.16--haef29d1_2'

    // input file can be VCF or BCF as is handled equally by bcftools
    input:
    tuple val(meta), file(bcf_allpos)

    output:
    tuple val(meta), path("${out_vcf}"),  emit: out_vcf

    script:
    out_vcf = "${meta.id}.vcf.gz"
    if (!params.report_ref_and_alt)
        """
        bcftools view -o ${out_vcf} \
            -O 'z' \
            -i 'GT="alt"' \
            '${bcf_allpos}'
        """
    else
        """
        bcftools view -o ${out_vcf} \
            -O 'z' \
            '${bcf_allpos}'
        """
}

process GATK_FINAL_VCF {
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'
    
    publishDir "${params.outdir}/${meta.id}/gatk/final_vcf", mode: 'copy', overwrite: true

    container 'quay.io/biocontainers/bcftools:1.16--haef29d1_2'

    // input file can be VCF or BCF as is handled equally by bcftools
    input:
    tuple val(meta), file(vcf_allpos)

    output:
    tuple val(meta), path("${out_vcf}"),  emit: out_vcf

    script:
    out_vcf = "${meta.id}.vcf.gz"
    if (!params.report_ref_and_alt)
        """
        bcftools view -o ${out_vcf} \
            -O 'z' \
            -i 'GT="alt"' \
            '${vcf_allpos}'
        """
    else
        """
        bcftools view -o ${out_vcf} \
            -O 'z' \
            '${vcf_allpos}'
        """
}
