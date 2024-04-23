process BCFTOOLS_MPILEUP {
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'

    // using package from conda-forge not bioconda (thus different from what underlies the biocontainers container) as there is a problem with lbgsl see https://github.com/samtools/bcftools/issues/1965
    conda 'conda-forge::gsl=2.7 bioconda::bcftools=1.17' 
    container 'quay.io/biocontainers/bcftools:1.17--h3cc50cf_1'

    input:
    tuple val(meta), path(sorted_reads_bam), path(sorted_reads_index), path(reference), path(reference_index)

    output:
    tuple val(meta), path("${mpileup_file}"),  emit: mpileup_file

    script:
    mpileup_file = "${meta.ID}.mpileup"
    """
    bcftools mpileup -o ${mpileup_file} \\
                     -O u \\
                     -f ${reference} \\
                     --min-BQ 20 \\
                     --annotate AD \\
                     ${sorted_reads_bam} 
    """
}

process BCFTOOLS_CALL {
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'

    // using package from conda-forge not bioconda (thus different from what underlies the biocontainers container) as there is a problem with lbgsl see https://github.com/samtools/bcftools/issues/1965
    conda 'conda-forge::gsl=2.7 bioconda::bcftools=1.17' 
    container 'quay.io/biocontainers/bcftools:1.17--h3cc50cf_1'

    input:
    tuple val(meta), file(mpileup_file)

    output:
    tuple val(meta), path("${vcf_allpos}"),  emit: vcf_allpos

    script:
    vcf_allpos = "${meta.ID}.vcf"
    """
    echo "${meta.ID}" > sample_name.txt
    bcftools call \\
        --output-type v \\
        -V indels \\
        -m \\
        --ploidy 1 \\
        ${mpileup_file} \\
        | bcftools reheader --samples sample_name.txt \\
        | bcftools view --output-file ${vcf_allpos} --output-type v
    """
}

process BCFTOOLS_FILTERING {
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'

    // using package from conda-forge not bioconda (thus different from what underlies the biocontainers container) as there is a problem with lbgsl see https://github.com/samtools/bcftools/issues/1965
    conda 'conda-forge::gsl=2.7 bioconda::bcftools=1.17' 
    container 'quay.io/biocontainers/bcftools:1.17--h3cc50cf_1'

    input:
    tuple val(meta), file(vcf_allpos)

    output:
    tuple val(meta), path("${filtered_vcf_allpos}"),  emit: filtered_vcf_allpos

    script:
    filtered_vcf_allpos = "${meta.ID}_filtered.vcf"
    """
    bcftools filter -o ${filtered_vcf_allpos} \\
                  ${params.VCF_filters} \\
                  --output-type v \\
                  ${vcf_allpos}
    """
}
process RAW_VCF {
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'
    
    publishDir "${params.outdir}/${meta.ID}/raw_vcf", mode: 'copy', overwrite: true

    // using package from conda-forge not bioconda (thus different from what underlies the biocontainers container) as there is a problem with lbgsl see https://github.com/samtools/bcftools/issues/1965
    conda 'conda-forge::gsl=2.7 bioconda::bcftools=1.17' 
    container 'quay.io/biocontainers/bcftools:1.17--h3cc50cf_1'

    // input file can be VCF or BCF as is handled equally by bcftools
    input:
    tuple val(meta), file(bcf_allpos)

    output:
    tuple val(meta), path("${out_vcf}"),  emit: out_vcf

    script:
    out_vcf = "${meta.ID}.vcf.gz"
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

    // using package from conda-forge not bioconda (thus different from what underlies the biocontainers container) as there is a problem with lbgsl see https://github.com/samtools/bcftools/issues/1965
    conda 'conda-forge::gsl=2.7 bioconda::bcftools=1.17' 
    container 'quay.io/biocontainers/bcftools:1.17--h3cc50cf_1'

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
