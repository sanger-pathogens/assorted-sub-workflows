process BASECALL {
    label 'cpu_1'
    label 'mem_4'
    label 'time_1'
    label 'gpu'

    tag "${params.barcode_kit_name}"

    container 'quay.io/sangerpathogens/cuda_dorado:0.9.1'
    
    input:
    path(pod5)

    output:
    path("calls.bam"), emit: called_channel

    script:
    def min_qscore = "${params.min_qscore == "" ? "" : "--min-qscore ${params.min_qscore}"}"

    def basecallCommand = "dorado basecaller ${params.basecall_model} --trim ${params.trim_adapters} ${min_qscore}"

    if (params.barcode_arrangement) {
        basecallCommand += " --kit-name ${params.barcode_kit_name} --barcode-arrangement ${params.barcode_arrangement} --barcode-sequences ${params.barcode_sequences} "
    } else {
        basecallCommand += " --kit-name ${params.barcode_kit_name} "
    }

    basecallCommand += "${pod5} > calls.bam"

    """
    ${basecallCommand}
    """
}

process DEMUX {
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    tag "${params.barcode_kit_name}"

    container 'quay.io/sangerpathogens/cuda_dorado:0.9.1'

    input:
    path(called_bam)

    output:
    path("barcodes/*.bam"), emit: called_channel

    script:
    """
    dorado demux --output-dir ./barcodes --no-classify ${called_bam}
    """
}

process DORADO_SUMMARY { 
    label 'cpu_1'
    label 'mem_500M'
    label 'time_1'
    
    tag "${params.barcode_kit_name}"

    publishDir path: "${params.outdir}/sequencing_summary/", mode: 'copy', overwrite: true, pattern: "summary.tsv"

    container 'quay.io/sangerpathogens/cuda_dorado:0.9.1'
    
    input:
    path(called_bam)

    output:
    path("summary.tsv"), emit: summary_channel

    script:
    """
    dorado summary ${called_bam} > summary.tsv
    """
}
