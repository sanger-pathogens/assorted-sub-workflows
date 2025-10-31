process BMTAGGER {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_10'
    label 'time_queue_from_normal'

    container 'quay.io/biocontainers/bmtagger:3.101--h470a237_4'

    input:
    tuple val(meta), path(first_read), path(second_read)

    output:
    tuple val(meta), path(first_read), path(second_read), emit: data_ch
    path(bmtagger_list), emit: bmtagger_list_ch

    script:
    bmtagger_list="${meta.ID}.bmtagger.list"
    """
    if [ "${params.run_trf}" == "false" && "${params.run_trimmomatic}" == "false" ]; then
        gunzip -c ${first_read} > ${meta.ID}_1.fastq
        gunzip -c ${second_read} > ${meta.ID}_2.fastq
        fq1=${meta.ID}_1.fastq
        fq2=${meta.ID}_2.fastq
    else
        fq1=${first_read}
        fq2=${second_read}
    fi

    # make tmp folder for bmtagger
    mkdir bmtagger_tmp
    # run bmtagger
    bmtagger.sh -b ${params.bmtagger_db}/${params.bmtagger_host}.bitmask -x ${params.bmtagger_db}/${params.bmtagger_host}.srprism -T bmtagger_tmp -q1 \\
	 -1 "\${fq1}" -2 "\${fq2}" \\
	 -o ${meta.ID}.bmtagger.list
    """
}