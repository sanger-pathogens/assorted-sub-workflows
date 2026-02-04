process COMBINE_BINS {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_100M'
    label 'time_30m'

    container 'quay.io/sangerpathogens/python-curl:3.11'

    input:
    tuple val(meta), path(bin)

    output:
    tuple val(meta), path(merged_assembly), emit: combined_bins

    script:
    merged_assembly = "${meta.ID}_merged_assembly.fa"
    """
    for i in \$(ls ${bin}); do cat ${bin}/\$i >> ${merged_assembly}; done
    """
}

process SPLIT_READS {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_12'

    container 'quay.io/biocontainers/pysam:0.23.3--py39hdd5828d_1'

    input:
    tuple val(meta), path(untrusted_contigs), path(mapped_reads)

    output:
    tuple val(meta), path("bin.*.{permissive,strict}_{1,2}.fastq.gz"),  emit: split_reads

    script:
    command = "${projectDir}/assorted-sub-workflows/mags_maker/reassemble/modules/bin/filter_reads_for_bin_reassembly.py"
    """
    ${command} ${untrusted_contigs} . ${params.strict_max} ${params.permissive_max} --mapped_reads ${mapped_reads}
    """
}

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

    script:
    command = "${projectDir}/assorted-sub-workflows/mags_maker/assemble/bin/rm_short_contigs.py"
    long_scaffolds = "long_${contigs.name}"
    min_contig_length = [params.maxbin2_min_contig, params.concoct_min_contig, params.metabat_min_contig].min()
    """
    ${command} ${min_contig_length} ${contigs} > ${long_scaffolds}
    """
}

process COLLECT_BINS {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_100M'
    label 'time_30m'

    container 'quay.io/sangerpathogens/python-curl:3.11'

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path(final_bin), emit: bin_ch

    script:
    final_bin = "${meta.ID}_reassembled_bin"
    """
    mkdir ${final_bin}
    cp ${contigs} ${final_bin}
    """
}

process RENAME_ORIGINAL {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_100M'
    label 'time_30m'

    container 'quay.io/sangerpathogens/python-curl:3.11'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path(final_name), emit: renamed_file

    script:
    final_name = "long_${meta.ID}_${fasta.baseName.replaceAll('\\.', '_')}_orgin.fasta"
    """
    cp ${fasta} ${final_name}
    """
}

process CHOOSE_BEST_BIN {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_100M'
    label 'time_30m'

    container 'quay.io/sangerpathogens/python-curl:3.11'

    input:
    tuple val(meta), path(bin), path(summary)

    output:
    tuple val(meta), path(final_bin), path("${meta.ID}_best_bins_summary.tsv"), emit: bins

    script:
    command = "${projectDir}/assorted-sub-workflows/mags_maker/reassemble/modules/bin/choose_best_bin.py"
    final_bin = "${meta.ID}_best_bins"
    """
    ${command} ${summary} ${bin} ${final_bin} --min-completeness ${params.min_completeness} --max-contamination ${params.max_contamination}
    """
}
