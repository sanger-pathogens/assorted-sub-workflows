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
    minimum_base_quality = "${params.minimum_base_quality == "default" ? "" : "--min-BQ ${params.minimum_base_quality}"}"
    """
    bcftools mpileup -o ${mpileup_file} \\
                     -O 'u' \\
                     ${minimum_base_quality} \\
                     -f ${reference} \\
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
    bcftools call --output ${vcf_allpos} \
                  --output-type 'v' \
                  --skip-variants indels \
                  --multiallelic-caller \
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
    bcftools filter --output-type 'u' \
                    --include 'GT!="0/1"' \
                    --soft-filter 'Het' \
                    '${vcf_allpos}' \
    | bcftools filter --output ${filtered_vcf_allpos} \
                      --output-type 'v' \
                      --include '${params.VCF_filters}' \
                      --soft-filter LowQual
    """
}

process BCFTOOLS_EXTRACT {
    label 'cpu_2'
    label 'mem_1'
    label 'time_1'

    publishDir "${params.outdir}/${meta.ID}/vcf/${filter_name}", mode: 'copy', overwrite: true

    // using package from conda-forge not bioconda (thus different from what underlies the biocontainers container) as there is a problem with lbgsl see https://github.com/samtools/bcftools/issues/1965
    conda 'conda-forge::gsl=2.7 bioconda::bcftools=1.17' 
    container 'quay.io/biocontainers/bcftools:1.17--h3cc50cf_1'

    input:
    tuple val(meta), file(vcf_allpos)
    val(filter)
    val(filter_name)

    output:
    tuple val(meta), path(filtered_vcf),  emit: filtered_vcf

    script:
    filtered_vcf = "${meta.ID}_${filter_name}.vcf.gz"
    """
    bcftools view --output ${filtered_vcf} \\
                  --output-type 'z' \\
                  --include '${filter}' \\
                  '${vcf_allpos}'
    """
}

process PUBLISH_VCF {
    label 'cpu_2'
    label 'mem_1'
    label 'time_30m'
    
    publishDir "${params.outdir}/${meta.ID}/vcf", mode: 'copy', overwrite: true

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
