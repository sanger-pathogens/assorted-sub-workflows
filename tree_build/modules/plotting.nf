process PLOT_TREE {
    label "cpu_1"
    label "mem_1"
    label "time_30m"

    publishDir "${params.outdir}/tree/", pattern: '*.png', mode: 'copy', overwrite: true

    container 'quay.io/sangerpathogens/rapidnj:2.3.2-c1'

    input:
    path(newick)

    output:
    path("*.png"), emit: plots

    script:
    """
    plot_tree.py ${newick} ${newick.baseName}.png
    """
}