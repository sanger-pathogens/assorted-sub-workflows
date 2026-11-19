process PASS_OR_FAIL_SYLPH {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_30m'

    input:
    tuple val(meta), path(report)

    output:
    tuple val(meta), env(pass_or_fail), emit: pass_or_fail

    script:
    """
    if [ "${params.sylph_profile}" = "true" ]; then
        top_genus_abun=\$(grep -P '\\|g__[^|]+\\t' ${report} | cut -f2 | sort -n | tail -1)
        top_species_abun=\$(grep -P '\\|s__[^|]+\\t' ${report} | cut -f2 | sort -n | tail -1)
        # no genus/species detected: treat as 0% so the sample fails
        top_genus_abun=\${top_genus_abun:-0}
        top_species_abun=\${top_species_abun:-0}
        genus_check=\$(echo "\$top_genus_abun < ${params.genus_abundance_threshold}" | bc)
        species_check=\$(echo "\$top_species_abun < ${params.species_abundance_threshold}" | bc)
        pass_or_fail=\$(if [ \$genus_check -eq 1 ] || [ \$species_check -eq 1 ]; then echo 'fail'; else echo 'pass'; fi)

    else
        pass_or_fail=\$(echo "NA")
    fi 
    """
}