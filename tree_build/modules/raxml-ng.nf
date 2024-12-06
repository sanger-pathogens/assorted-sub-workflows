process RAXML_NG {
    label 'cpu_2'
    label 'mem_16'
    label 'time_12'

    conda 'bioconda::raxml-ng=1.1.0'
    container 'quay.io/biocontainers/raxml-ng:1.1.0--h22e3c99_1'

    publishDir "${params.outdir}/tree", mode: 'copy', overwrite: true, pattern: "*.support"

    input:
    tuple path(msa), path(constant_sites_freq)

    output:
    path("*.support"), emit: tree

    script:
    """
    constant_sites_freq=\$(cat ${constant_sites_freq})
    raxml_model="${params.base_model}+ASC_STAM{\${constant_sites_freq}}"
    raxml-ng --check --msa ${msa} --model \${raxml_model}
    raxml-ng --parse --msa ${msa} --model \${raxml_model}
    raxml-ng --all \\
        --msa ${msa}.raxml.rba \\
        --model \${raxml_model} \\
        --tree ${params.tree_search} \\
        --bs-trees ${params.bootstrap_trees} \\
        --prefix ${msa.baseName} \\
        --threads auto{${task.cpus}} \\
        ${params.raxmlng_args}
    """
}
