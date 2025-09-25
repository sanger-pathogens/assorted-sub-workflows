// define min contig length
min_contig_length = [params.maxbin2_min_contig, params.concoct_min_contig, params.metabat_min_contig].min()

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
    ${command} ${min_contig_length} ${contigs} > ${long_scaffolds} 2> >(grep "Warning:" > remove_small_contigs.err)
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
    ${command} ${min_contig_length} ${contigs} > ${long_scaffolds}
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
    ${command} *scaffolds --min_contig ${min_contig_length} > ${final_contigs}
    """
}
