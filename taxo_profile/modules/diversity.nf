process CALC_DIVERSITY {
    // tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_12'


    // publishDir "${params.outdir}/${meta.ID}/sylph/", pattern: "*.sylsp", mode: 'copy', overwrite: true, enabled: params.save_sylph_sketches

    // container 'quay.io/biocontainers/sylph:0.8.1--ha6fb395_0'

    input:
    path(sylph_profiles)

    output:
    path("diversity.tsv"), emit: diversity

    script:
    def diversity = "${moduleDir}/../scripts/calculate_diverity.py"
    """
    ${diversity} -i *_sylph_profile.tsv \
        --outdir ./ \
        --taxonomic-abundance-threshold ${params.taxonomic_abundance_threshold} 
    """
}