process GENERATE_STATS {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'
    
    container 'quay.io/sangerpathogens/metawrap_qc_python:1.0'

    input:
    tuple val(meta), path("trimmed_read_1.fq"), path("trimmed_read_2.fq"), path("clean_read_1.fq.gz"), path("clean_read_2.fq.gz"), path("host_read_1.fq.gz"), path("host_read_2.fq.gz"), path("original_read_1.fq.gz"), path("original_read_2.fq.gz")

    output:
    path(stats_file), emit: stats_ch

    script:
    stats_file="${meta.ID}_stats.csv"
    """
    # get read numbers (bash quickest way)
    trimmed_1_reads=\$((`cat trimmed_read_1.fq | wc -l` / 4))
    trimmed_2_reads=\$((`cat trimmed_read_2.fq | wc -l` / 4))
    clean_1_reads=\$((`zcat clean_read_1.fq.gz | wc -l` / 4))
    clean_2_reads=\$((`zcat clean_read_2.fq.gz | wc -l` / 4))
    host_1_reads=\$((`zcat host_read_1.fq.gz | wc -l` / 4))
    host_2_reads=\$((`zcat host_read_2.fq.gz | wc -l` / 4))
    original_1_reads=\$((`zcat original_read_1.fq.gz | wc -l` / 4))
    original_2_reads=\$((`zcat original_read_2.fq.gz | wc -l` / 4))
    trimmed_reads_total=\$((\${trimmed_1_reads} + \${trimmed_2_reads}))
    clean_reads_total=\$((\${clean_1_reads} + \${clean_2_reads}))
    host_reads_total=\$((\${host_1_reads} + \${host_2_reads}))
    original_reads_total=\$((\${original_1_reads} + \${original_2_reads}))

    # generate stats
    generate_stats.py --sample-id ${meta.ID} --host-reads \${host_reads_total} --non-host-reads \${clean_reads_total} --total-trimmed-reads \${trimmed_reads_total} --total-original-reads \${original_reads_total} > ${meta.ID}_stats.csv
    """
}
