process GTDBTK {
    tag "${meta.ID}"

    def mem_label = params.temp_file_storage ? 'mem_16' : 'mem_120'
    
    label "cpu_8"
    label "time_12"
    label mem_label

    container 'quay.io/biocontainers/gtdbtk:2.4.1--pyhdfd78af_1'

    publishDir mode: 'copy', path: "${params.outdir}/gtdbtk/"

    input:
    tuple val(meta), path(fastas, stageAs: "fastas/*")

    output:
    tuple val(meta), path(report_tsv), emit: results

    script:
    def scratch_dir_base = params.temp_file_storage
    def use_scratch = (scratch_dir_base != null)

    def scratch_setup = use_scratch ? """
        [[ -d "${scratch_dir_base}" ]] || { echo "Creating scratch base: ${scratch_dir_base}"; mkdir -p "${scratch_dir_base}"; }
        SCRATCH_DIR=\$(mktemp -d -p "${scratch_dir_base}" gtdbtk_temp_XXXXXXXX)
        echo "GTDB-Tk scratch dir: \$SCRATCH_DIR" >&2
        cleanup() { rm -rf "\$SCRATCH_DIR"; }
        trap cleanup EXIT

    """ : """
        echo "GTDB-Tk running WITHOUT --scratch_dir" >&2
        SCRATCH_DIR=""
    """

    report_tsv = "${meta.ID}_gtdbtk_summary.tsv"

    """
    set -euo pipefail
    export GTDBTK_DATA_PATH="${params.gtdbtk_db}"

    ${scratch_setup}

    scratch_flag=""
    if [ -n "\$SCRATCH_DIR" ]; then
      scratch_flag="--scratch_dir \$SCRATCH_DIR"
    fi

    gtdbtk classify_wf --genome_dir fastas -x ${params.fasta_ext} --skip_ani_screen --cpus ${task.cpus} --out_dir gtdbtk_outdir \$scratch_flag

    cp gtdbtk_outdir/gtdbtk.bac*.summary.tsv ${report_tsv}
    """

}

