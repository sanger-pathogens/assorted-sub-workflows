process BCFTOOLS_MPILEUP {
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'

    container 'quay.io/biocontainers/bcftools:1.16--haef29d1_2'

    input:
    tuple val(meta), path(sorted_reads), path(reference), path(reference_index)

    output:
    tuple val(meta), path("${mpileup_file}"),  emit: mpileup_file

    script:
    mpileup_file = "${meta.id}.mpileup"
    """
    bcftools mpileup -o ${mpileup_file} \
                     -O 'u' \
                     -f ${reference} \
                     ${sorted_reads} 
    """
}

process BCFTOOLS_CALL {
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'

    container 'quay.io/biocontainers/bcftools:1.16--haef29d1_2'

    input:
    tuple val(meta), file(mpileup_file)

    output:
    tuple val(meta), path("${vcf_allpos}"),  emit: vcf_allpos

    script:
    vcf_allpos = "${meta.id}.vcf"
    """
    bcftools call -o ${vcf_allpos} \
        -O 'v' \
        -V indels \
        -m \
        '${mpileup_file}'
    """
}

process BCFTOOLS_FILTERING {
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
                  -i '${params.VCF_filters}' \
                  '${vcf_allpos}'
    """
}
process RAW_VCF {
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'
    
    publishDir "${params.outdir}/${meta.id}/raw_vcf", mode: 'copy', overwrite: true

    container 'quay.io/biocontainers/bcftools:1.16--haef29d1_2'

    // input file can be VCF or BCF as is handled equally by bcftools
    input:
    tuple val(meta), file(bcf_allpos)

    output:
    tuple val(meta), path("${out_vcf}"),  emit: out_vcf

    script:
    out_vcf = "${meta.id}.vcf.gz"
    if (params.only_report_alts)
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

process FINAL_VCF {
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'
    
    publishDir "${params.outdir}/${meta.ID}/final_vcf", mode: 'copy', overwrite: true

    container 'quay.io/biocontainers/bcftools:1.16--haef29d1_2'

    // input file can be VCF or BCF as is handled equally by bcftools
    input:
    tuple val(meta), file(vcf_allpos)

    output:
    tuple val(meta), path("${out_vcf}"),  emit: out_vcf

    script:
    out_vcf = "${meta.ID}.vcf.gz"
    if (params.only_report_alts)
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
