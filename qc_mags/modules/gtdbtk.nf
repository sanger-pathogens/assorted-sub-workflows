process GTDBTK {
    tag "${meta.ID}"
    label((params.fast_taxonomy) ? (params.batch_size > 100 ? "cpu_4" : "cpu_2") : (params.batch_size > 100 ? "cpu_64" : "cpu_8"))
    label((params.fast_taxonomy || params.temp_file_storage) ? "mem_64" : "mem_120")
    label "time_queue_from_normal"

    container  'quay.io/biocontainers/gtdbtk:2.8.0--pyhdfd78af_0'

    publishDir mode: 'copy', path: "${params.outdir}/${qc_stage}/gtdbtk/", enabled: !(params.skip_raw_reports)

    input:
    tuple val(meta), path(fastas, stageAs: "fastas/*")
    val(qc_stage)

    output:
    tuple val(meta), path(report_tsv), emit: results

    script:
    report_tsv = "${meta.ID}_gtdbtk_summary.tsv"
    place_flag = params.fast_taxonomy ? "" : "--place_species"

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

    gtdbtk classify_wf --genome_dir fastas -x ${params.fasta_ext} ${place_flag} --cpus ${task.cpus} --out_dir gtdbtk_outdir \$scratch_flag

    awk 'FNR==1 && NR!=1 {next} {print}' gtdbtk_outdir/gtdbtk.bac*.summary.tsv \$(ls gtdbtk_outdir/gtdbtk.ar*.summary.tsv 2>/dev/null) > ${report_tsv}
    """

}

