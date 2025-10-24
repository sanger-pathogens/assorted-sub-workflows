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
    gtdbtk classify_wf --genome_dir fastas -x ${params.fasta_ext} --skip_ani_screen --cpus ${task.cpus} --out_dir gtdbtk_outdir

    cp gtdbtk_outdir/gtdbtk.bac*.summary.tsv ${report_tsv}
    """
}

process skip_GTDBTK {
    tag "${meta.ID}"
    label "cpu_1"
    label "mem_1"
    label "time_1"

    container 'quay.io/sangerpathogens/pandas:2.2.1'
    
    publishDir mode: 'copy', path: "${params.outdir}/gtdbtk/"

    input:
    tuple val(meta), path(fastas, stageAs: "fastas/*")

    output:
    tuple val(meta), path(report_tsv), emit: results

    script:
    skeleton_tsv = "${projectDir}/assorted-sub-workflows/qc_isolates/assets/gtdbtk_silence.tsv"
    report_tsv = "${meta.ID}_gtdbtk_summary.tsv"
    """
    #!/usr/bin/env python
    import os
    import pandas as pd

    inputs = os.listdir(fastas)
    skeleton = pd.read_csv("${skeleton_tsv}", sep='\\t')
    headers = list(skeleton.columns)

    rows = []
    for input in inputs:
        name = os.path.splitext(input)[0]
        row = [name] + ['na']*len(headers)
        rows.append(row)

    df = pd.DataFrame(rows, columns=headers)
    df.to_csv("${report_tsv}", sep='\\t', index=False)
    """
}