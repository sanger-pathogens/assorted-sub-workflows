process FILTER_HOST_READS {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'
    
    container 'quay.io/sangerpathogens/metawrap_qc_python:1.0'

    input:
    tuple val(meta), path(first_read), path(second_read), path(bmtagger_list)

    output:
    tuple val(meta), path(first_read), path(second_read), emit: data_ch
    tuple val(meta), path(clean_1), path(clean_2), emit: cleaned_ch

    script:
    clean_1="${meta.ID}_clean_1.fastq"
    clean_2="${meta.ID}_clean_2.fastq"
    """
    # filter out host reads
    filter_reads.py -b ${bmtagger_list} -r ${first_read} --skip-human-reads > ${clean_1}
    filter_reads.py -b ${bmtagger_list} -r ${second_read} --skip-human-reads > ${clean_2}
    """
}

process FILTER_HOST_READS_SINGLE {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'
    
    container 'quay.io/sangerpathogens/metawrap_qc_python:1.0'

    input:
    tuple val(meta), path(fastq), path(bmtagger_list)

    output:
    tuple val(meta), path(fastq), emit: data_ch
    tuple val(meta), path(clean_fastq), emit: cleaned_ch

    script:
    clean_fastq = "${meta.ID}_unpaired_clean.fastq"
    """
    filter_reads.py -b ${bmtagger_list} -r ${fastq} --skip-human-reads > ${clean_fastq}
    """
}

process GET_HOST_READS {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_4'  //TODO: Need 5GB?
    label 'time_12'

    if (params.publish_host_reads) { publishDir path: "${params.outdir}/host_reads", mode: 'copy', overwrite: true, pattern: "*_host*.fastq.gz" }
    container 'quay.io/sangerpathogens/metawrap_qc_python:1.0'

    input:
    tuple val(meta), path(first_read), path(second_read), path(bmtagger_list)

    output:
    tuple val(meta), path(host_1), path(host_2), emit: host_ch

    script:
    host_1="${meta.ID}_host_1.fastq.gz"
    host_2="${meta.ID}_host_2.fastq.gz"
    """
    # filter out host reads
    filter_reads.py -b ${bmtagger_list} -r ${first_read} --get-human-reads > ${meta.ID}_host_1.fastq
    filter_reads.py -b ${bmtagger_list} -r ${second_read} --get-human-reads > ${meta.ID}_host_2.fastq
    # compress
    pigz ${meta.ID}_host_1.fastq ${meta.ID}_host_2.fastq
    """
}

process GET_HOST_READS_SINGLE {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_4'  //TODO: Need 5GB?
    label 'time_12' 

    if (params.publish_host_reads) { publishDir path: "${params.outdir}/host_reads", mode: 'copy', overwrite: true, pattern: "*_host*.fastq.gz" }
    container 'quay.io/sangerpathogens/metawrap_qc_python:1.0'

    input:
    tuple val(meta), path(fastq), path(bmtagger_list)

    output:
    tuple val(meta), path(host_fastq), emit: host_ch

    script:
    host_fastq = "${meta.ID}_unpaired_host.fastq.gz"
    """
    filter_reads.py -b ${bmtagger_list} -r ${fastq} --get-human-reads > ${meta.ID}_unpaired_host.fastq
    pigz ${meta.ID}_unpaired_host.fastq
    """
}
