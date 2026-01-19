process BASECALL {
    label 'cpu_1'
    label 'mem_4'
    label 'time_1'
    label 'gpu'

    tag "${params.barcode_kit_name}"

    container 'quay.io/sangerpathogens/cuda_dorado:1.3.1'
    
    input:
    path(pod5)

    output:
    path("calls.bam"), emit: called_channel

    script:
    def barcodeArgs = ""
    if (params.barcode_arrangement) {
        barcodeArgs = "--kit-name ${params.barcode_kit_name} " +
                    "--barcode-arrangement ${params.barcode_arrangement} " +
                    "--barcode-sequences ${params.barcode_sequences} "
    } else if (params.barcode_kit_name) {
        barcodeArgs = "--kit-name ${params.barcode_kit_name} "
    }

    def min_qscore_args = params.min_qscore == "" ? "" : "--min-qscore ${params.min_qscore}"

    def methylation_models_args = params.modified_bases_models ? "--modified-bases-models  ${params.modified_bases_models}" : ""
    

    def basecallCommand =
        "dorado basecaller ${params.model} --trim ${params.trim_adapters} " +
        "${min_qscore_args} ${barcodeArgs} ${methylation_models_args} ${pod5} > calls.bam"
    
    """
    ${basecallCommand}
    """
}

process DEMUX {
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    tag "${params.barcode_kit_name}"

    container 'quay.io/sangerpathogens/cuda_dorado:1.3.1'

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

    container 'quay.io/sangerpathogens/cuda_dorado:1.3.1'
    
    input:
    path(called_bam)

    output:
    path("summary.tsv"), emit: summary_channel

    script:
    """
    dorado summary ${called_bam} > summary.tsv
    """
}
