process SYLPH_SKETCH {
    tag "${meta.ID}"
    label 'cpu_4'
    label 'mem_4'
    label 'time_queue_from_small'

    publishDir "${params.outdir}/${meta.ID}/sylph/", pattern: "*.sylsp", mode: 'copy', overwrite: true, enabled: params.save_sylph_sketches

    container 'quay.io/biocontainers/sylph:1.0.0--hb42e459_0'

    input:
    tuple val(meta), path(read_1), path(read_2)

    output:
    tuple val(meta), path("${meta.ID}.paired.sylsp"), emit: sketch

    script:
    """
    sylph sketch -t ${task.cpus} -1 ${read_1} -2 ${read_2} -k ${params.sylph_k} -S ${meta.ID} -d ./
    """
}

process SYLPH_PROFILE {
    tag "batch_${batch_index}"
    label 'cpu_8'
    label 'mem_32'
    label 'time_queue_from_normal'

    publishDir params.outdir, saveAs: { filename -> "${filename - ~/_sylph_profile\.tsv$/}/sylph/${filename}" }, pattern: "*_sylph_profile.tsv", mode: 'copy', overwrite: true

    container 'quay.io/biocontainers/sylph:1.0.0--hb42e459_0'

    input:
    tuple val(batch_index), path(sketches)

    output:
    tuple val(batch_index), path("batch_${batch_index}.tsv"), emit: batch_report
    path("*_sylph_profile.tsv"), optional: true, emit: sylph_reports

    script:
    """
    sylph profile -t ${task.cpus} -k ${params.sylph_k} ${params.sylph_profile_options} \\
        -o batch_${batch_index}.tsv ${sketches} ${params.sylph_db}

    # rows are grouped by sample (first column, the sketch name given by -S)
    awk -F'\\t' 'NR==1 {header=\$0; next}
                \$1 != sample {if (out) close(out); sample=\$1; out=sample"_sylph_profile.tsv"; print header > out}
                {print > out}' batch_${batch_index}.tsv
    """
}

process SYLPHTAX_TAXPROF {
    tag "batch_${batch_index}"
    label 'cpu_1'
    label 'mem_4'
    label 'time_queue_from_small'

    publishDir params.outdir, saveAs: { filename -> "${file(filename).baseName}/sylph/${filename}" }, pattern: "*.sylphmpa", mode: 'copy', overwrite: true

    container 'quay.io/biocontainers/sylph-tax:1.9.2--pyhdfd78af_0'

    input:
    tuple val(batch_index), path(sylph_report), path(sylph_tax_metadata)

    output:
    path("*.sylphmpa"), emit: sylphtax_mpa_reports

    script:
    """
    sylph-tax taxprof "${sylph_report}" -t "${sylph_tax_metadata}"
    """
}
