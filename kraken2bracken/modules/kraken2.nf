process KRAKEN2 {
    tag "${meta.ID}"
    label 'cpu_4'
    label 'mem_8'
    label 'time_12'

    publishDir "${params.outdir}/${meta.ID}/kraken2", mode: 'copy', overwrite: true, pattern: "*.tsv"

    container 'quay.io/biocontainers/kraken2:2.1.3--pl5321hdcf5f25_0'

    input:
    tuple val(meta), path(read_1), path(read_2), path(kraken2_db)

    output:
    tuple val(meta), path("*kraken_report.tsv"),  emit: kraken2_report
    tuple val(meta), path("*kraken_sample_report.tsv"),  emit: kraken2_sample_report

    script:
    kraken2_id = "${meta.ID}".replaceAll('#','_')
    memory_mapping = (params.memory_mapping ? "--memory-mapping" : "")
    """
    kraken2 --db "${kraken2_db}" \
            --threads ${task.cpus} \
            --output "${meta.ID}_kraken_report.tsv" \
            --use-names \
            --report "${meta.ID}_kraken_sample_report.tsv" \
            --report-zero-counts \
            --report-minimizer-data \
            --paired "${read_1}" "${read_2}" \
            ${memory_mapping}
    """
}

process KRAKEN2_GET_CLASSIFIED {
    tag "${meta.ID}"
    label 'cpu_4'
    label 'mem_8'
    label 'time_12'

    publishDir "${params.outdir}/${meta.ID}/kraken2", mode: 'copy', overwrite: true, pattern: "*.tsv"
    publishDir "${params.outdir}/${meta.ID}/kraken2", mode: 'copy', overwrite: true, pattern: "classified.fastq"

    container 'quay.io/biocontainers/kraken2:2.1.3--pl5321hdcf5f25_0'

    input:
    tuple val(meta), path(read_1), path(read_2), path(kraken2_db)

    output:
    tuple val(meta), path("*kraken_report.tsv"),  emit: kraken2_report
    tuple val(meta), path("*kraken_sample_report.tsv"),  emit: kraken2_sample_report
    tuple val(meta), path("*_classified.fastq"),  emit: classified_reads
    tuple val(meta), path("*_unclassified.fastq"),  emit: unclassified_reads

    script:
    kraken2_id = "${meta.ID}".replaceAll('#','_')
    memory_mapping = (params.memory_mapping ? "--memory-mapping" : "")
    """
    if [ "${params.memory_mapping}" = true ]; then
        echo "du -sh \$PWD"
    fi
    kraken2 --db "${kraken2_db}" \
            --threads ${task.cpus} \
            --classified-out "${kraken2_id}#_classified.fastq" --unclassified-out "${kraken2_id}#_unclassified.fastq" \
            --output "${meta.ID}_kraken_report.tsv" \
            --use-names \
            --report "${meta.ID}_kraken_sample_report.tsv" \
            --report-zero-counts \
            --report-minimizer-data \
            --paired "${read_1}" "${read_2}" \
            ${memory_mapping}        
    if [ "${params.memory_mapping}" = true ]; then
        echo "du -sh \$PWD"
    fi    
    """
}

process COMPRESS_READS {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    publishDir "${params.outdir}/${meta.ID}/kraken2", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_classified.fastq.gz"),  emit: classified_reads
    tuple val(meta), path("*_unclassified.fastq.gz"),  emit: unclassified_reads

    script:
    """
    gzip -fk *.fastq
    """
}
