process BOWTIE2_HOST_REMOVAL {
    tag "${meta.ID}"
    label 'cpu_8'
    label 'mem_10'
    label 'time_queue_from_long'

    container 'ghcr.io/wtsi-npg/p4:0.48.5'

    input:
    tuple val(meta), path(read_1), path(read_2)

    output:
    tuple val(meta), path("${meta.ID}_hostfree_1.fastq.gz"), path("${meta.ID}_hostfree_2.fastq.gz"), emit: reads
    path("${meta.ID}_host_removal_stats.tsv"), emit: stats
    path("${meta.ID}.bowtie2.log"), emit: log

    script:
    """
    # interleave the pair and normalise read names to <name>/1, <name>/2
    paste <(zcat -f ${read_1} | paste - - - -) <(zcat -f ${read_2} | paste - - - -) \\
    | awk -F'\\t' -v OFS='\\n' '{split(\$1,a," "); sub(/\\/1\$/,"",a[1]); print a[1]"/1",\$2,"+",\$4,a[1]"/2",\$6,"+",\$8}' \\
    | bowtie2 -p ${task.cpus} --local --very-sensitive-local --reorder -x ${params.host_bowtie2_index} -U - 2> ${meta.ID}.bowtie2.log \\
    | perl -nale 'BEGIN{\$,="\\t"} if(/^\\@/){print; next} if(\$F[0]=~s{/([12])\$}{}){ \$F[1] = \$1==1 ? \$F[1]|0x41 : \$F[1]|0x81 } print @F' \\
    | samtools fixmate -u - - \\
    | samtools fastq -f 12 -N -0 /dev/null -s /dev/null -c 4 -@ 2 \\
        -1 ${meta.ID}_hostfree_1.fastq.gz -2 ${meta.ID}_hostfree_2.fastq.gz -
    cat ${meta.ID}.bowtie2.log >&2

    reads_in=\$(head -1 ${meta.ID}.bowtie2.log | awk '{print \$1}')
    aligned=\$(grep -E 'aligned (exactly 1 time|>1 times)' ${meta.ID}.bowtie2.log | awk '{s+=\$1} END{print s+0}')
    kept=\$(( \$(zcat ${meta.ID}_hostfree_1.fastq.gz | wc -l) / 4 * 2 ))
    printf "Sample_id\\treads_in\\treads_aligned_human\\treads_kept\\tpct_removed\\n%s\\t%s\\t%s\\t%s\\t%s\\n" ${meta.ID} \$reads_in \$aligned \$kept \\
        \$(awk -v i=\$reads_in -v k=\$kept 'BEGIN{printf "%.4f", i ? 100*(i-k)/i : 0}') > ${meta.ID}_host_removal_stats.tsv
    """
}
