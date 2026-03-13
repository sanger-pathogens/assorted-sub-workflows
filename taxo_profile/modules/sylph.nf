process SYLPH_SKETCH {
    tag "${meta.ID}"
    label 'cpu_2'
    label 'mem_2'
    label 'time_from_queue_small'

    publishDir "${params.outdir}/${meta.ID}/sylph/", pattern: "*.sylsp", mode: 'copy', overwrite: true, enabled: params.save_sylph_sketches

    container 'gitlab-registry.internal.sanger.ac.uk/sanger-pathogens/docker-images/sylph:0.8.1--ha6fb395_0'
    errorStrategy 'terminate'

    input:
    tuple val(meta), path(read_1), path(read_2)

    output:
    tuple val(meta), path("${meta.ID}.paired.sylsp"), emit: sketch

    script:
    """
    sylph sketch -t ${task.cpus} -1 ${read_1} -2 ${read_2} -k ${params.sketch_size} -S ${meta.ID} -d ./
    """
}

process SYLPH_PROFILE {
    tag "${meta.ID}"
    label 'cpu_2'
    label 'mem_20'
    label 'time_from_queue_small'

    publishDir "${params.outdir}/${meta.ID}/sylph/", pattern: "*.tsv", mode: 'copy', overwrite: true

    container 'gitlab-registry.internal.sanger.ac.uk/sanger-pathogens/docker-images/sylph:0.8.1--ha6fb395_0'
    errorStrategy 'terminate'

    input:
    tuple val(meta), path(sketch)

    output:
    tuple val(meta), path("${meta.ID}_sylph_profile.tsv"), emit: sylph_report

    script:
    def estimate_unknown = params.sylph_estimate_unknown ? (params.sylph_read_seq_id ? "-u --read-seq-id ${params.sylph_read_seq_id}" : "-u") : ""
    """
    sylph profile -t ${task.cpus} -o ${meta.ID}_sylph_profile.tsv -k ${params.sketch_size} ${sketch} ${params.sylph_db} ${estimate_unknown}
    """
}

process SYLPH_QUERY {
    tag "${meta.ID}"
    label 'cpu_2'
    label 'mem_20'
    label 'time_from_queue_small'

    publishDir "${params.outdir}/${meta.ID}/sylph/", pattern: "*.tsv", mode: 'copy', overwrite: true

    container 'gitlab-registry.internal.sanger.ac.uk/sanger-pathogens/docker-images/sylph:0.8.1--ha6fb395_0'
    errorStrategy 'terminate'

    input:
    tuple val(meta), path(sketch)

    output:
    tuple val(meta), path("${meta.ID}_sylph_profile.tsv"), emit: sylph_report

    script:
    """
    sylph query ${sketch} ${params.sylph_db} -t ${task.cpus} -o ${meta.ID}_sylph_profile.tsv
    """
}

process SYLPHTAX_TAXPROF {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_2'
    label 'time_from_queue_small'

    publishDir "${params.outdir}/${meta.ID}/sylph/", pattern: "*.sylphmpa", mode: 'copy', overwrite: true

    container 'quay.io/biocontainers/sylph-tax:1.7.0--pyhdfd78af_0'
    errorStrategy 'terminate'

    input:
    tuple val(meta), path(sylph_report)
    path sylph_tax_metadata

    output:
    tuple val(meta), path("${meta.ID}_sylphtax_profile.sylphmpa") , emit: sylphtax_mpa_report

    script:
    """
    metadata_file=\$(basename "${sylph_tax_metadata}")
    normalized_report="${meta.ID}_sylph_tax_input.tsv"
    python3 - <<'PYTHON'
import csv
from pathlib import Path

input_path = Path("${sylph_report}")
output_path = Path("${meta.ID}_sylph_tax_input.tsv")
columns = [
    "Sample_file",
    "Genome_file",
    "Taxonomic_abundance",
    "Sequence_abundance",
    "Adjusted_ANI",
    "Eff_cov",
    "ANI_5-95_percentile",
    "Eff_lambda",
    "Lambda_5-95_percentile",
    "Median_cov",
    "Mean_cov_geq1",
    "Containment_ind",
    "Naive_ANI",
    "kmers_reassigned",
    "Contig_name",
]

with input_path.open(newline="") as src, output_path.open("w", newline="") as dst:
    reader = csv.DictReader(src, delimiter="	")
    writer = csv.DictWriter(dst, fieldnames=columns, delimiter="	", extrasaction="ignore")
    writer.writeheader()
    for row in reader:
        writer.writerow({col: row.get(col, "") for col in columns})
PYTHON
    sylph-tax taxprof "${normalized_report}" -t "\${metadata_file}"
    mv ${meta.ID}.sylphmpa ${meta.ID}_sylphtax_profile.sylphmpa
    """
}


process SYLPH_SUMMARIZE {
    label 'cpu_1'
    label 'mem_4'
    label 'time_queue_from_small'

    container 'quay.io/sangerpathogens/pandas:2.2.1'
    errorStrategy 'terminate'

    publishDir "${params.outdir}/sylph/", mode: 'copy', overwrite: true

    input:
    path(sylph_reports)

    output:
    path("references.txt"), emit: references
    path("sylph_summary.tsv"), emit: sylph_summary

    script:
    // Filter once with thresholds.
    """
    ${workflow.projectDir}/assorted-sub-workflows/sylph_refset/bin/sylph_summarize.py \
        --reports ${sylph_reports} \
        --ani ${params.sylph_ani} \
        --cov ${params.sylph_cov} \
        --ani-column Naive_ANI \
        --cov-column Eff_cov \
        --out-references references.txt \
        --out-summary sylph_summary.tsv
    """
}
