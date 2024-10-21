process PASS_OR_FAIL_FASTQC {
    label 'cpu_1'
    label 'mem_1'
    label 'time_30m'

    container 'ubuntu:24.04'

    input:
    tuple val(meta), path(read_1_zip), path(read_2_zip)

    output:
    tuple val(meta), val(pass_or_fail), emit: pass_or_fail

    script:
    """
    unzip ${read_1_zip}
    unzip ${read_2_zip}
    pass_or_fail=\$(if grep 'FAIL' */summary.txt; then echo 'fail'; else echo 'pass'; fi)
    """
}

process PASS_OR_FAIL_K2B {
    label 'cpu_1'
    label 'mem_1'
    label 'time_30m'

    container 'ubuntu:24.04'

    input:
    tuple val(meta), path(report)

    output:
    tuple val(meta), val(pass_or_fail), emit: pass_or_fail

    script:
    """
    top_genus_abun=\$(grep -P 'G\t' ${report} | cut -f1 | sort -n | tail -1)
    top_species_abun=\$(grep -P 'S\t' ${report} | cut -f1 | sort -n | tail -1)
    if [ top_genus_abun -lt ${params.genus_abundance_threshold} ] || [ top_species_abun -lt ${params.species_abundance_threshold} ]
    then
        pass_or_fail = 'fail'
    else
        pass_or_fail = 'pass'
    fi
    """
}
