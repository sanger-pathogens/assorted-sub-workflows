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
    label 'time_1'

    container 'quay.io/sangerpathogens/python-curl:3.11'

    input:
    tuple val(meta), path(untrusted_contigs), path(sam)

    output:
    tuple val(meta), path("bin.*.{permissive,strict}_{1,2}.fastq.gz"),  emit: split_reads

    script:
    command = "${projectDir}/modules/reassemble/bin/filter_reads_for_bin_reassembly.py"
    """
    ${command} ${untrusted_contigs} . ${params.strict_max} ${params.permissive_max} --sam ${sam}
    """
}
