process CALC_DIVERSITY {
    label 'cpu_1'
    label 'mem_1'
    label 'time_12'

    publishDir "${params.outdir}/", pattern: "diversity_estimates.tsv", mode: 'copy', overwrite: true
    container 'quay.io/sangerpathogens/python_graphics:1.1.7'

    input:
    path(sylph_profiles)

    output:
    path("diversity_estimates.tsv"), emit: diversity

    script:
    def diversity = "${moduleDir}/../scripts/calculate_diverity.py"
    """
    ${diversity} -i *_sylph_profile.tsv \
        --outdir ./ \
        --taxonomic-abundance-threshold ${params.taxonomic_abundance_threshold} 
    """
}