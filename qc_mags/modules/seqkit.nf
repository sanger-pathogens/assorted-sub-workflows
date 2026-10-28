include { fastaExt; genomeName } from './helper_functions.nf'

process SEQKIT {
    tag "${meta.ID}"
    label "cpu_1"
    label "mem_1"
    label "time_30m"

    container  'quay.io/biocontainers/seqkit:2.10.0--h9ee0642_0'

    publishDir mode: 'copy', path: "${params.outdir}/seqkit/"

    input:
    tuple val(group_key), val(meta), path(fasta)

    output:
    tuple val(group_key), val(meta), path(finalName), emit: results

    script:
    // MDMCLEANER output is <name>_filtered_kept_contigs; isolates arrive uncleaned, so both get the same output name
    genome_name = genomeName(fasta).replaceAll(/_filtered_kept_contigs$/, '')
    finalName = "cleaned_${genome_name}_filtered_kept_contigs.${fastaExt()}"
    """
    seqkit seq ${fasta} -m ${params.min_contig} > ${finalName}

    if [[ ! -s ${finalName} ]]; then
        echo "Assembly (fasta) ${finalName} is empty after filtering contigs < ${params.min_contig} bp" 1>&2 && exit 3
    fi
    """
}