process COLLATE_STATS_BMTAGGER {
    tag "all samples"
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    publishDir "${params.outdir}/preprocessing_summary_stats", mode: 'copy', overwrite: true, pattern: "*_statistics.csv"
    
    input:
    path(stats_files)

    output:
    path(read_removal_stats_file), emit: host_reads_stats_ch

    script:
    read_removal_stats_file="read_removal_statistics.csv"
    """
    echo "Sample_id,Total_host_reads,Total_non_host_reads,Total_trimmed_reads,host_reads_%,non_host_reads_%,Total_original_reads,reads_trimmed_%" > "${read_removal_stats_file}"
    cat *_stats.csv >> "${read_removal_stats_file}"
    """
}

process COLLATE_STATS_TRIMMOMATIC {
    tag "all samples"
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    publishDir "${params.outdir}/preprocessing_summary_stats", mode: 'copy', overwrite: true, pattern: "*_statistics.csv"
    
    input:
    path(stats_files)

    output:
    path(trimming_stats_file), emit: trimming_stats_ch

    script:
    trimming_stats_file="trimmomatic_statistics.csv"
    """
    echo "Sample_id,Input_read_pairs,Both_surviving_reads,Both_surviving_read_percent,Forward_only_surviving_reads,Forward_only_surviving_read_percent,Reverse_only_surviving_reads,Reverse_only_surviving_read_percent,Dropped_reads,Dropped_read_percent" > "${trimming_stats_file}"

    for file in *_trimmomatic_summary.csv; do

        sample_id=\$(basename "\$file" | sed 's/_trimmomatic_summary.*//')
        input_read_pairs=\$(grep "Input Read Pairs" "\$file" | awk '{print \$4}' | xargs)
        both_surviving_reads=\$(grep "Both Surviving Reads" "\$file" | awk '{print \$4}' | xargs)
        both_surviving_read_percent=\$(grep "Both Surviving Read Percent" "\$file" | awk '{print \$5}' | xargs)
        forward_only_surviving_reads=\$(grep "Forward Only Surviving Reads" "\$file" | awk '{print \$5}' | xargs)
        forward_only_surviving_read_percent=\$(grep "Forward Only Surviving Read Percent" "\$file" | awk '{print \$6}' | xargs)
        reverse_only_surviving_reads=\$(grep "Reverse Only Surviving Reads" "\$file" | awk '{print \$5}' | xargs)
        reverse_only_surviving_read_percent=\$(grep "Reverse Only Surviving Read Percent" "\$file" | awk '{print \$6}' | xargs)
        dropped_reads=\$(grep "Dropped Reads" "\$file" | awk '{print \$3}' | xargs)
        dropped_read_percent=\$(grep "Dropped Read Percent" "\$file" | awk '{print \$4}' | xargs)

        echo "\${sample_id},\${input_read_pairs},\${both_surviving_reads},\${both_surviving_read_percent},\${forward_only_surviving_reads},\${forward_only_surviving_read_percent},\${reverse_only_surviving_reads},\${reverse_only_surviving_read_percent},\${dropped_reads},\${dropped_read_percent}" >> "${trimming_stats_file}"
    
    done
    """
}