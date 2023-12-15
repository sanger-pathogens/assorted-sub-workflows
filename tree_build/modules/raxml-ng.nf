process BUILD_TREE {
    label 'cpu_2'
    label 'mem_16'
    label 'time_12'

    container '/software/pathogen/images/raxml-ng-1.1.0.simg'

    publishDir "${params.outdir}/tree", mode: 'copy', overwrite: true, pattern: "*.support"

    input:
    tuple path(msa), val(conscount)

    output:
    path("*.support"), emit: tree_channel

    script:
    conssites = conscount.text.trim()
    raxml_model="${params.raxml_base_model}+ASC_STAM{${conssites}}"
    """
    raxml-ng --check --msa ${msa} --model ${raxml_model}
    raxml-ng --parse --msa ${msa} --model ${raxml_model}
    raxml-ng --all --msa ${msa}.raxml.rba --model ${raxml_model} \
    --tree pars{10} --bs-trees 200 \
    --threads ${params.raxml_threads}
    """
}
