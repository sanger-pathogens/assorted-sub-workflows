process BINETTE {
    tag "${meta.ID}"
    label 'cpu_8'
    label 'mem_64'
    label 'time_12'

    container 'quay.io/biocontainers/binette:1.2.1--pyh106432d_1'

    publishDir mode: 'copy', path: "${params.outdir}/${meta.ID}/binette/", pattern: "${report_txt}"

    input:
    tuple val(meta), path(comebin_bins, stageAs: 'comebin_bins'),
                     path(semibin2_bins, stageAs: 'semibin2_bins'),
                     path(metacat_bins, stageAs: 'metacat_bins'),
                     path(maxbin2_bins, stageAs: 'maxbin2_bins'),
                     path(assembly)

    output:
    tuple val(meta), path("final_bins/*"), path(report_txt), emit: results

    script:
    report_txt = "${meta.ID}_binette_quality_report.tsv"
    """
    binette \\
        --bin_dirs comebin_bins semibin2_bins metacat_bins maxbin2_bins \\
        --fasta_extensions .fasta .fa .fna \\
        --contigs ${assembly} \\
        --checkm2_db ${params.checkm2_db} \\
        -o binette_out -t ${task.cpus}

    mv binette_out/final_bins final_bins
    mv binette_out/final_bins_quality_reports.tsv ${report_txt}
    """
}
