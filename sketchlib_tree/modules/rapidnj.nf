process RAPIDNJ {
    label "cpu_4"
    label "mem_8"
    label "time_1"

    container 'quay.io/sangerpathogens/rapidnj:2.3.2-c1'

    publishDir "${params.outdir}/tree/", pattern: '*.nwk', mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(total_tsv)

    output:
    tuple val(meta), path("*.nwk"), emit: tree

    script:
    ani_tree_tools = "${projectDir}/assorted-sub-workflows/sketchlib_tree/bin/ani_tree_tools.py"
    """
    ${ani_tree_tools} --dist_tsv_path ${total_tsv} --meta_ID ${meta.ID} --core_accession --build_tree
    """
}