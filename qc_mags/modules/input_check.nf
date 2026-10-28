process INPUT_FASTA_CHECK {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    container 'quay.io/sangerpathogens/pandas:2.2.1'

    input:
    tuple val(meta), path(fastas, stageAs: 'input/*')

    output:
    tuple val(meta), path('checked/*'),             emit: fastas
    path(fixes),                                    emit: fixes

    script:
    command = "${projectDir}/assorted-sub-workflows/qc_mags/bin/check_input_fasta.py"
    fixes = "${meta.ID}_input_fasta_fixes.tsv"
    """
    ${command} \\
        --batch ${meta.ID} \\
        --fixes ${fixes} \\
        --outdir checked \\
        input/*
    """
}
