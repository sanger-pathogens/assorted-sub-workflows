process BCFTOOLS_MPILEUP {
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'

    // using package from conda-forge not bioconda (thus different from what underlies the biocontainers container) as there is a problem with lbgsl see https://github.com/samtools/bcftools/issues/1965
    conda 'conda-forge::gsl=2.7 bioconda::bcftools=1.17' 
    container 'quay.io/biocontainers/bcftools:1.17--h3cc50cf_1'

    input:
    tuple val(meta), path(sorted_reads), path(reference), path(reference_index)

    output:
    tuple val(meta), path("${mpileup_file}"),  emit: mpileup_file

    script:
    mpileup_file = "${meta.ID}.mpileup"
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

    // using package from conda-forge not bioconda (thus different from what underlies the biocontainers container) as there is a problem with lbgsl see https://github.com/samtools/bcftools/issues/1965
    conda 'conda-forge::gsl=2.7 bioconda::bcftools=1.17' 
    container 'quay.io/biocontainers/bcftools:1.17--h3cc50cf_1'

    input:
    tuple val(meta), file(mpileup_file)

    output:
    tuple val(meta), path("${vcf_allpos}"),  emit: vcf_allpos

    script:
    vcf_allpos = "${meta.id}_raw.vcf"
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
    bcftools view -o ${filtered_vcf_allpos} \
                  -O 'v' \
                  -i '${params.bcftools_vcf_filter}' \
                  '${vcf_allpos}'
    """
}

process BCFTOOLS_VIEW {
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'
    
    publishDir "${params.outdir}/${meta.id}/bcftools/raw_vcf", mode: 'copy', overwrite: true, pattern: "*_raw.vcf.gz"
    publishDir "${params.outdir}/${meta.id}/bcftools/final_vcf", mode: 'copy', overwrite: true, pattern: "*_filtered.vcf.gz"

    // using package from conda-forge not bioconda (thus different from what underlies the biocontainers container) as there is a problem with lbgsl see https://github.com/samtools/bcftools/issues/1965
    conda 'conda-forge::gsl=2.7 bioconda::bcftools=1.17' 
    container 'quay.io/biocontainers/bcftools:1.17--h3cc50cf_1'

    // input file can be VCF or BCF as is handled equally by bcftools
    input:
    tuple val(meta), file(bcf_allpos)

    output:
    tuple val(meta), path("${out_vcf}"),  emit: out_vcf

    script:
    vcf_prefix = bcf_allpos.simpleName
    out_vcf = "${vcf_prefix}.vcf.gz"
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
