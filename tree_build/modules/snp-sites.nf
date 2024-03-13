process SNP_SITES{
    label 'cpu_2'
    label 'mem_100M'
    label 'time_12'

    conda 'bioconda::snp-sites=2.5.1'
    container 'quay.io/biocontainers/snp-sites-2.5.1--hed695b0_0'

    publishDir "${params.outdir}/snp_aln", mode: 'copy', overwrite: true, pattern: "*.snp.aln"

    input:
    path(msa)

    output:
    tuple path(output_snpaln), path(output_conscount), emit: snp_aln_channel

    script:
    output_snpaln="${msa}.snp.aln"
    output_conscount="${msa}.conscount"
    """
    snp-sites -o ${output_snpaln} ${msa}
    snp-sites -C ${msa} | tr ',' '/' > ${msa}.conscount
    """

}