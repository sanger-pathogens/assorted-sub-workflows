process BOWTIE_INDEX {
    label 'cpu_4'
    label 'mem_32'
    label 'time_queue_from_normal'

    container 'quay.io/sangerpathogens/bowtie2-samtools:1.1-c1'

    input:
    tuple val(meta), path(subset_fasta)

    output:
    tuple val(meta), path("*.bt2*"), emit: bowtie_index

    script:
    """
    bowtie2-build $subset_fasta ${meta.id}_gtdb_subset.bt2 --threads 4 --large-index
    """
}

process BOWTIE2SAMTOOLS {
    tag "${meta.id}"
    label 'cpu_4'
    label 'mem_4'
    label 'time_queue_from_normal'

    container 'quay.io/sangerpathogens/bowtie2-samtools:1.1-c1'

    if (params.bowtie2_samtools_only_abundance_estimation) { publishDir path: "${params.results_dir}", mode: 'copy', overwrite: true, pattern: "*.sorted.bam" }
    input:
    tuple val(meta), path(first_read), path(second_read), val(btidx)
    val threads

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), emit: bam_file
    tuple val(meta), path(first_read), path(second_read), emit: trimmed_fastqs
    path("${meta.id}_mapping_rate.csv"), emit: map_rate_ch

    script:
    """
    bowtie2 -p $threads -x ${btidx} -1 $first_read -2 $second_read | samtools sort -@ $threads -o ${meta.id}".sorted.bam"
    mapping_rate=\$(grep "overall alignment rate" .command.err | awk '{ print \$1 }')
    echo ${meta.id},\${mapping_rate} > ${meta.id}_mapping_rate.csv
    """
}

process GET_OVERALL_MAPPING_RATE {
    
    publishDir "${params.results_dir}/mapping_rates/", mode: 'copy', overwrite: true, pattern: 'mapping_rates.csv', saveAs: { filename -> "${workflow.start}_mapping_rates.csv" }
    
    input:
    path(mapping_rate)

    output:
    path("mapping_rates.csv")

    script:
    """
    echo "Sample,Mapping_Rate" > mapping_rates.csv
    cat *_mapping_rate.csv >> mapping_rates.csv
    """
}