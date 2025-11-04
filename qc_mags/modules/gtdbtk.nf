process GTDBTK {
    tag "${meta.ID}"
    label "cpu_8"
    label "time_12"
    label params.temp_file_storage ? "mem_16" : "mem_120"


    container  'quay.io/biocontainers/gtdbtk:2.4.1--pyhdfd78af_1'

    publishDir mode: 'copy', path: "${params.outdir}/${qc_stage}/gtdbtk/", enabled: !(params.skip_raw_reports)

    input:
    tuple val(meta), path(fastas, stageAs: "fastas/*")
    val(qc_stage)

    output:
    tuple val(meta), path(report_tsv), emit: results

    script:
    report_tsv = "${meta.ID}_gtdbtk_summary.tsv"

    """
    if [[ -n ${params.temp_file_storage} && "${params.temp_file_storage}" != "null" ]]; then
        [[ -d "${params.temp_file_storage}" ]] || { echo "Creating scratch base: ${params.temp_file_storage}"; mkdir -p "${params.temp_file_storage}"; }
        SCRATCH_DIR=\$(mktemp -d -p "${params.temp_file_storage}" gtdbtk_temp_XXXXXXXX)
        echo "GTDB-Tk scratch dir: \$SCRATCH_DIR" >&2
        cleanup() { rm -rf "\$SCRATCH_DIR"; }
        trap cleanup EXIT
    else
        echo "GTDB-Tk running WITHOUT --scratch_dir" >&2
        SCRATCH_DIR=""
    fi

    scratch_flag=""
    if [ -n "\$SCRATCH_DIR" ]; then
      scratch_flag="--scratch_dir \$SCRATCH_DIR"
    fi

    export GTDBTK_DATA_PATH="${params.gtdbtk_db}"

    gtdbtk classify_wf --genome_dir fastas -x ${params.fasta_ext} --skip_ani_screen --cpus ${task.cpus} --out_dir gtdbtk_outdir \$scratch_flag

    cp gtdbtk_outdir/gtdbtk.bac*.summary.tsv ${report_tsv}

    """

}

