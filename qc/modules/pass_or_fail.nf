process PASS_OR_FAIL_FASTQC {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_30m'

    input:
    tuple val(meta), path(read_1_zip), path(read_2_zip)

    output:
    tuple val(meta), env(pass_or_fail), emit: pass_or_fail

    script:
    """
    unzip ${read_1_zip}
    unzip ${read_2_zip}
    pass_or_fail=\$(if grep 'FAIL' */summary.txt; then echo 'fail'; else echo 'pass'; fi)
    """
}

process PASS_OR_FAIL_K2B {
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
    top_genus_abun=\$(grep -P 'G\t' ${report} | cut -f1 | sort -n | tail -1)
    top_species_abun=\$(grep -P 'S\t' ${report} | cut -f1 | sort -n | tail -1)
    genus_check=\$(echo "\$top_genus_abun < ${params.genus_abundance_threshold}" | bc)
    species_check=\$(echo "\$top_species_abun < ${params.species_abundance_threshold}" | bc)
    pass_or_fail=\$(if [ \$genus_check -eq 1 ] || [ \$species_check -eq 1 ]; then echo 'fail'; else echo 'pass'; fi)
    """
}
