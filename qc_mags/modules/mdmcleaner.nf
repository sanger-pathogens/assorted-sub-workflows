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
    // MDMcleaner names its outputs after the genome name
    genome_name = genomeName(fasta)
    """
    # remove colons common in JGI/IMG style genome names and name them something neutral to prefer barnapp error
    mkdir -p header_safe
    if gzip -t ${fasta} 2>/dev/null; then
        zcat ${fasta} | sed '/^>/ s/:/__COLON__/g' | gzip > header_safe/qcgenome.${fasta_ext}.gz
        staged=header_safe/qcgenome.${fasta_ext}.gz
    else
        sed '/^>/ s/:/__COLON__/g' ${fasta} > header_safe/qcgenome.${fasta_ext}
        staged=header_safe/qcgenome.${fasta_ext}
    fi

    mdmcleaner set_configs --db_basedir ${params.mdmcleaner_db} --threads ${task.cpus}
    mdmcleaner clean -i \$staged -o mdmcleaner_output -t ${task.cpus} -c mdmcleaner.config --fast_run

    # restore the genome's own name on MDMcleaner's output folder and files
    if [ -d mdmcleaner_output/qcgenome ]; then
        for file in mdmcleaner_output/qcgenome/qcgenome_*; do
            mv "\$file" "mdmcleaner_output/qcgenome/${genome_name}_\${file##*/qcgenome_}"
        done
        mv mdmcleaner_output/qcgenome mdmcleaner_output/${genome_name}
    fi

    for file in mdmcleaner_output/*/*.fasta.gz; do
        zcat "\$file" | sed '/^>/ s/__COLON__/:/g' | gzip > "\$file.tmp"
        mv "\$file.tmp" "\$file"
    done

    # Rename files to match desired extension as mdmcleaner is always fasta
    for file in mdmcleaner_output/*/*.fasta.gz; do
        new_name=\$(echo "\$file" | sed "s|\\.fasta\\.gz\$|.${output_fasta_ext}|")
        mv "\$file" "\$new_name"
    done

    mv mdmcleaner_output/*/*_filtered_kept_contigs.${output_fasta_ext} .
    """
}
