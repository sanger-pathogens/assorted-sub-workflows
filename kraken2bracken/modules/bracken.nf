process SEARCH_LIBRARY {
    label 'cpu_4'
    label 'mem_8'
    label 'time_12'

    container 'quay.io/biocontainers/bracken:2.8--py310h0dbaff4_1'

    input:
    path(kraken2_db)

    output:
    tuple path(kraken2_db), path("database.kraken"),  emit: kraken2_library

    script:
    """
    kraken2 --db=${kraken2_db} \
            --threads=${task.cpus} \
            <( find -L ${kraken2_db}/library \\( -name "*.fna" -o -name "*.fasta" -o -name "*.fa" \\) -exec cat {} + ) \
            > database.kraken
    """
}

process KMER2READ_DIST {
    label 'cpu_4'
    label 'mem_8'
    label 'time_12'

    container 'quay.io/biocontainers/bracken:2.8--py310h0dbaff4_1'

    input:
    tuple path(kraken2_db), path(kraken2_library)

    output:
    path("database${params.read_len}mers.kraken"),  emit: kmer2read_distr

    script:
    """
    kmer2read_distr --seqid2taxid ${kraken2_db}/seqid2taxid.map \
                    --taxonomy ${kraken2_db}/taxonomy \
                    --kraken ${kraken2_library} \
                    --output database${params.read_len}mers.kraken \
                    -k ${params.kmer_len} \
                    -l ${params.read_len} \
                    -t ${task.cpus}
    """
}

process GENERATE_KMER_DISTRIBUTION {
    label 'cpu_4'
    label 'mem_8'
    label 'time_12'

    publishDir "${params.outdir}/bracken_build", mode: 'copy', overwrite: true, pattern: "*.kmer_distrib"

    container 'quay.io/biocontainers/bracken:2.8--py310h0dbaff4_1'

    input:
    path(kmer2read_distr)

    output:
    path("*.kmer_distrib"),  emit: kmer_distrib

    script:
    """
    generate_kmer_distribution.py -i database${params.read_len}mers.kraken \
                                  -o database${params.read_len}mers.kmer_distrib
    """
}

process BRACKEN {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    publishDir "${params.outdir}/${meta.ID}/bracken", mode: 'copy', overwrite: true

    container 'quay.io/biocontainers/bracken:2.8--py310h0dbaff4_1'

    input:
    tuple val(meta), path(kraken2_sample_report), path(kmer_distrib)

    output:
    tuple val(meta), path("*.bracken"),  emit: bracken_report
    tuple val(meta), path("*_kraken_sample_report_bracken_*.tsv"),  emit: kraken_style_bracken_report

    script:
    // biocontainer doesn't expose the script est_abundance.py on PATH :(
    """
    /usr/local/bin/src/est_abundance.py -i ${kraken2_sample_report} \
                                        -k ${kmer_distrib} \
                                        -o ${meta.ID}.bracken \
                                        -l ${params.bracken_classification_level} \
                                        -t ${params.threshold}
    """
}
