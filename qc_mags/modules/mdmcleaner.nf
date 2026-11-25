include { fastaExt; genomeName } from './helper_functions.nf'

process MDMCLEANER {
    tag "${meta.ID}"
    label "cpu_16"
    label "mem_8"
    label "time_12"

    container  'quay.io/biocontainers/mdmcleaner:0.8.7--pyh7cba7a3_0'

    publishDir enabled: params.publish_mdmcleaner , mode: 'copy', path: "${params.outdir}/mdmcleaner/"

    input:
    tuple val(group_key), val(meta), path(fasta)

    output:
    tuple val(group_key), val(meta), path("*.${output_fasta_ext}"), emit: results

    script:
    fasta_ext = fastaExt()
    output_fasta_ext = "${fasta_ext}.gz"
    // MDMcleaner runs on a copy called qcgenome; the kept contigs are written back out under the genome name
    genome_name = genomeName(fasta)
    """
    # MDMcleaner's barrnap step fails on any file name containing "fail"/"error", and MDMcleaner splits contig
    # names on ':' (JGI/IMG style). Run it on a neutrally named copy with ':' swapped out, then restore both.
    zcat -f ${fasta} | sed '/^>/ s/:/__COLON__/g' > qcgenome.${fasta_ext}

    mdmcleaner set_configs --db_basedir ${params.mdmcleaner_db} --threads ${task.cpus}
    mdmcleaner clean -i qcgenome.${fasta_ext} -o mdmcleaner_output -t ${task.cpus} -c mdmcleaner.config --fast_run

    zcat mdmcleaner_output/qcgenome/qcgenome_filtered_kept_contigs.fasta.gz \
        | sed '/^>/ s/__COLON__/:/g' \
        | gzip > ${genome_name}_filtered_kept_contigs.${output_fasta_ext}
    """
}
