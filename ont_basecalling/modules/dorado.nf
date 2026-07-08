process BASECALL {
    tag "${meta.ID}"

    label 'cpu_1'
    label 'mem_2'
    label 'time_12'
    label 'gpu'

    container 'quay.io/sangerpathogens/cuda_dorado:1.3.1'
    
    input:
    tuple val(meta), path(pod5)

    output:
    tuple val(meta), path("calls.bam"), emit: called_channel
    tuple val(meta), path(pod5), path("${meta.ID}_basecalling_complete.txt"), emit: basecalling_complete_channel

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

    touch ${meta.ID}_basecalling_complete.txt
    """
}

process BASECALL_LEGACY {
    tag "${meta.ID}"

    label 'cpu_1'
    label 'mem_2'
    label 'time_12'
    label 'gpu'

    container 'quay.io/sangerpathogens/cuda_dorado:0.9.6'
    
    input:
    tuple val(meta), path(pod5)

    output:
    tuple val(meta), path("calls.bam"), emit: called_channel
    tuple val(meta), path(pod5), path("${meta.ID}_basecalling_complete.txt"), emit: basecalling_complete_channel

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

    touch ${meta.ID}_basecalling_complete.txt
    """
}

process DEMUX {
    tag "${meta.ID}"

    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    tag "${params.barcode_kit_name}"

    container 'quay.io/sangerpathogens/cuda_dorado:1.3.1'

    input:
    tuple val(meta), path(called_bam)

    output:
    tuple val(meta), path("barcodes/*.bam"), emit: called_channel

    script:
    """
    dorado demux --output-dir ./barcodes --no-classify ${called_bam}
    """
}

process DORADO_SUMMARY {
    tag "${meta.ID}"
    
    label 'cpu_1'
    label 'mem_500M'
    label 'time_1'
    
    tag "${params.barcode_kit_name}"

    publishDir path: "${params.outdir}/qc/", mode: 'copy', overwrite: true

    container 'quay.io/sangerpathogens/cuda_dorado:1.3.1'
    
    input:
    tuple val(meta), path(called_bam)

    output:
    tuple val(meta), path(summary), emit: summary_channel

    script:
    summary = "${meta.ID}_sequencing_summary.txt"
    """
    dorado summary ${called_bam} > ${summary}
    """
}
