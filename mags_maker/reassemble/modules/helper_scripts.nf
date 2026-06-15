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

    container 'quay.io/sangerpathogens/python-curl:3.11'

    input:
    tuple val(meta), path(untrusted_contigs), path(sam)

    output:
    tuple val(meta), path("bin.*.{permissive,strict}_{1,2}.fastq.gz"),  emit: split_reads

    script:
    command = "${projectDir}/assorted-sub-workflows/mags_maker/reassemble/modules/bin/filter_reads_for_bin_reassembly.py"
    """
    ${command} ${untrusted_contigs} . ${params.strict_max} ${params.permissive_max} --sam ${sam}
    """
}

process COLLECT_BINS {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_100M'
    label 'time_30m'

    container 'quay.io/sangerpathogens/python-curl:3.11'

    input:
    tuple val(meta), path(contigs), val(rename_map)

    output:
    tuple val(meta), path(final_bin), emit: bin_ch

    script:
    final_bin = "${meta.ID}_reassembled_bin"
    rename_cmds = rename_map
        .collect { staged, final_name -> "[ -f '${staged}' ] && mv '${staged}' '${final_name}'" }
        .join('\n')
    """
    mkdir ${final_bin}
    ${rename_cmds}
    cp ${contigs} ${final_bin}/
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
