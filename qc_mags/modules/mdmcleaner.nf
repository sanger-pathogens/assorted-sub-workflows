process MDMCLEANER {
    tag "${meta.ID}"
    label "cpu_32"
    label "mem_8"
    label "time_12"

    container  'quay.io/biocontainers/mdmcleaner:0.8.7--pyh7cba7a3_0'

    publishDir mode: 'copy', path: "${params.outdir}/mdmcleaner/"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.${params.fasta_ext}.gz"), emit: results

    script:
    """
    mdmcleaner set_configs --db_basedir ${params.mdmcleaner_db} --threads ${task.cpus}
    mdmcleaner clean -i ${fasta}/*${params.fasta_ext} -o mdmcleaner_output -t ${task.cpus} -c mdmcleaner.config --fast_run


    # Rename files to match desired extension as mdmcleaner is always fasta
    for file in mdmcleaner_output/*/*_filtered_kept_contigs.fasta.gz; do
        new_name=\$(echo "\$file" | sed "s|\\.fasta\\.gz\$|\\.${params.fasta_ext}.gz|")
        mv "\$file" "\$new_name"
    done

    mv mdmcleaner_output/*/*_filtered_kept_contigs.${params.fasta_ext}.gz .
    """
}