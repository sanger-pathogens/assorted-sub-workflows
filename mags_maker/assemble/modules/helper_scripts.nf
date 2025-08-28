process REMOVE_SMALL_CONTIGS {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_100M'
    label 'time_30m'

    container 'quay.io/sangerpathogens/python-curl:3.11'

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path(long_scaffolds), emit: long_contigs
    path('remove_small_contigs.err'), emit: warning_log

    script:
    command = "${projectDir}/assorted-sub-workflows/mags_maker/assemble/bin/rm_short_contigs.py"
    long_scaffolds = "${meta.ID}_long.scaffolds"
    """
    ${command} ${params.min_contig} ${contigs} > ${long_scaffolds} 2> >(grep "Warning:" > remove_small_contigs.err)
    """
}

process FIX_MEGAHIT_CONTIG_NAMING {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_100M'
    label 'time_30m'

    container 'quay.io/sangerpathogens/python-curl:3.11'

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path(long_scaffolds), emit: long_contigs

    script:
    command = "${projectDir}/assorted-sub-workflows/mags_maker/assemble/bin/fix_megahit_contig_naming.py"
    long_scaffolds = "${meta.ID}_long.scaffolds"
    """
    ${command} ${params.min_contig} ${contigs} > ${long_scaffolds}
    """
}

process SORT_CONTIGS {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_100M'
    label 'time_30m'

    container 'quay.io/sangerpathogens/python-curl:3.11'

    input:
    tuple val(meta), path("${meta.ID}_long?.scaffolds")

    output:
    tuple val(meta), path(final_contigs), emit: sorted_contigs

    script:
    command = "${projectDir}/assorted-sub-workflows/mags_maker/assemble/bin/sort_contigs.py"
    final_contigs = "${meta.ID}.contigs"
    """
    ${command} *scaffolds --min_contig ${params.min_contig} > ${final_contigs}
    """
}
