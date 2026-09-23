process REMOVE_SMALL_CONTIGS {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_100M'
    label 'time_12'

    container 'quay.io/sangerpathogens/python-curl:3.11'

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path(long_scaffolds), emit: long_contigs
    path('remove_small_contigs.err'), emit: warning_log

    script:
    command = "${projectDir}/assorted-sub-workflows/mags_maker/assemble/bin/rm_short_contigs.py"
    long_scaffolds = "${meta.ID}_long.scaffolds"
    min_contig_length = [params.maxbin2_min_contig, params.concoct_min_contig, params.metabat_min_contig].min()
    """
    ${command} ${min_contig_length} ${contigs} > ${long_scaffolds} 2> >(grep "Warning:" > remove_small_contigs.err)
    """
}

process FIX_MEGAHIT_CONTIG_NAMING {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_100M'
    label 'time_12'

    container 'quay.io/sangerpathogens/python-curl:3.11'

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path(long_scaffolds), emit: long_contigs

    script:
    command = "${projectDir}/assorted-sub-workflows/mags_maker/assemble/bin/fix_megahit_contig_naming.py"
    long_scaffolds = "${meta.ID}_long.scaffolds"
    min_contig_length = [params.maxbin2_min_contig, params.concoct_min_contig, params.metabat_min_contig].min()
    """
    ${command} ${min_contig_length} ${contigs} > ${long_scaffolds}
    """
}

process SORT_CONTIGS {
    tag "${meta.ID}"
    label 'cpu_1'
    // mem_100M was sized on a 0.55 GB/mate dev sample (~25 MB assembly, 44 MB
    // peak RSS). sort_contigs.py holds the WHOLE assembly in memory twice -
    // parse_fasta() builds a dict of every contig, then combine_fastas() builds
    // a second of those passing the length filter - so peak is ~2x the assembly
    // in Python strings plus per-object overhead. On a 968 MB assembly that is
    // GBs: SRR25448209 died here with exit 130 on all three escalations
    // (100 -> 200 -> 400 MB) and was then dropped silently by errorStrategy
    // 'ignore'. Five samples were lost this way before it was found.
    label 'mem_8'
    label 'time_12'

    publishDir mode: 'copy', path: "${params.outdir}/${meta.ID}/pre_binning", saveAs: { filename -> "${meta.ID}_unbinned_contigs.fasta" }, enabled: params.publish_unbinned_contigs

    container 'quay.io/sangerpathogens/python-curl:3.11'

    input:
    tuple val(meta), path("${meta.ID}_long?.scaffolds")

    output:
    tuple val(meta), path(final_contigs), emit: sorted_contigs

    script:
    command = "${projectDir}/assorted-sub-workflows/mags_maker/assemble/bin/sort_contigs.py"
    final_contigs = "${meta.ID}.contigs"
    min_contig_length = [params.maxbin2_min_contig, params.concoct_min_contig, params.metabat_min_contig].min()
    """
    ${command} *scaffolds --min_contig ${min_contig_length} > ${final_contigs}
    """
}
