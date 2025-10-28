process SEQTK_MERGEPE {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    container 'quay.io/biocontainers/seqtk:1.4--he4a0461_2'

    input:
    tuple val(meta), path(read_1), path(read_2)

    output:
    tuple val(meta), path("${meta.ID}_interleaved.fq")

    script:

    """
    seqtk mergepe ${read_1} ${read_2} > ${meta.ID}_interleaved.fq
    """
}

process SEQTK_SPLIT{
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    container 'quay.io/biocontainers/seqtk:1.4--he4a0461_2'

    publishDir enabled: params.debug_preproc_output, "${params.outdir}/${meta.ID}/preprocessing/", mode: "copy", pattern:"*.fq"

    input:
    tuple val(meta), path(interleaved_fq)
    val(suffix)

    output:
    tuple val(meta), path("${meta.ID}_${suffix}_split_1.fq"), path("${meta.ID}_${suffix}_split_2.fq")

    script:
    
    """
    ## this produces .fa files as output but the content of the files are proper .fq
    ## cmd is: seqtk -n [NUM OF FILES] [OUTPUT PREFIX] path/to/input

    seqtk split -n 2 ${meta.ID} ./${interleaved_fq}

    mv ${meta.ID}.00001.fa ${meta.ID}_${suffix}_split_1.fq
    mv ${meta.ID}.00002.fa ${meta.ID}_${suffix}_split_2.fq
    """
}