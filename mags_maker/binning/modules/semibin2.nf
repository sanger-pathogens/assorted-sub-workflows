process SEMIBIN2 {
    tag "${meta.ID}"
    label 'cpu_4'
    label 'mem_4'
    label 'time_12'

    container 'quay.io/biocontainers/semibin:2.4.1--pyhdfd78af_0'

    input:
    tuple val(meta), path(bam), path(bai), path(assembly)

    output:
    tuple val(meta), path(bins_out), emit: bins

    script:
    bins_out = "semibin2_out/output_bins_unzipped"
    """
    SemiBin2 single_easy_bin -i ${assembly} -b ${bam} -o semibin2_out \\
        -t ${task.cpus} --engine cpu --environment ${params.semibin2_environment}

    mkdir -p ${bins_out}
    for f in semibin2_out/output_bins/*.fa.gz; do
        gunzip -c "\$f" > "${bins_out}/\$(basename "\$f" .gz)"
    done
    """
}
