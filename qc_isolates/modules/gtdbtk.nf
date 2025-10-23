process GTDBTK {
    tag "${meta.ID}"
    label "cpu_8"
    label "mem_120"
    label "time_12"

    container  'quay.io/biocontainers/gtdbtk:2.4.1--pyhdfd78af_1'

    publishDir mode: 'copy', path: "${params.outdir}/gtdbtk/"

    input:
    tuple val(meta), path(fastas, stageAs: "fastas/*")

    output:
    tuple val(meta), path(report_tsv), emit: results

    script:
    report_tsv = "${meta.ID}_gtdbtk_summary.tsv"
    """
    export GTDBTK_DATA_PATH="${params.gtdbtk_db}"

    gtdbtk classify_wf --genome_dir fastas -x ${params.fasta_ext} --skip_ani_screen --cpus ${task.cpus} --out_dir gtdbtk_outdir --scratch_dir ${tmp_dir}

    cp gtdbtk_outdir/gtdbtk.bac*.summary.tsv ${report_tsv}
    """

}