process MDMCLEANER {
    tag "${meta.ID}"
    label "cpu_32"
    label "mem_8"
    label "time_12"

    container  'quay.io/biocontainers/mdmcleaner:0.8.7--pyh7cba7a3_0'

    publishDir mode: 'copy', path: "${params.outdir}/mdmcleaner/"

    input:
    tuple val(meta), path(fastas, stageAs: "fastas/*")

    output:
    tuple val(meta), path("*.${output_fasta_ext}"), emit: results

    script:
    fasta_ext = params.fasta_ext.replaceAll(/^\./, '')
    output_fasta_ext = "${fasta_ext}.gz"
    """
    mdmcleaner set_configs --db_basedir ${params.mdmcleaner_db} --threads ${task.cpus}
    mdmcleaner clean -i fastas/*.${fasta_ext} -o mdmcleaner_output -t ${task.cpus} -c mdmcleaner.config --fast_run

    # Rename files to match desired extension as mdmcleaner is always fasta
    for file in mdmcleaner_output/*/*.fasta.gz; do
        new_name=\$(echo "\$file" | sed "s|\\.fasta\\.gz\$|.${output_fasta_ext}|")
        mv "\$file" "\$new_name"
    done

    mv mdmcleaner_output/*/*_filtered_kept_contigs.${output_fasta_ext} .
    """
}