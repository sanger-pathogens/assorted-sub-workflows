include { CONVERT_FAST5_TO_POD5; MERGE_POD5 } from './modules/pod5.nf'
include { BASECALL; DEMUX; DORADO_SUMMARY   } from './modules/dorado.nf'
include { PYCOQC                            } from './modules/pycoqc.nf'
include { CONVERT_TO_FASTQ; PUBLISH_BAMS    } from './modules/samtools.nf'

def validateSingleFormat(listOfFormats){
    if (listOfFormats.size() != 1) {
        log.error("Multiple signal filetypes ${listOfFormats} found in '${params.raw_read_dir}'. Please separate filetypes into distinct directories and process independently.")
    }
}

def validateBarcodeParams() {
    def kitName = params.barcode_kit_name != null && params.barcode_kit_name.trim()
    def barcodeArgs = params.barcode_arrangement != null && params.barcode_sequences != null
    
    if (barcodeArgs && !kitName) {
        log.error("--barcode_kit_name must also be specified this is the name in the .toml.")
    }

    if (kitName && !barcodeArgs || (kitName && barcodeArgs)) {
        return // Valid configurations: either barcode_kit_name alone OR all three together
    }
    //if above return doesn't run throw exeception
    log.error("Must specify either only --barcode_kit_name or all three: --barcode_kit_name, --barcode_arrangement, and --barcode_sequences.")
}


workflow ONT_BASECALLING{  
    take:
    raw_read_signal_files

    main:
    validateBarcodeParams()

    raw_read_signal_files.map{ raw_read_signal_file ->
        tuple(raw_read_signal_file.extension, raw_read_signal_file)
    }
    | groupTuple
    | map { format, files -> 
            validateSingleFormat([format])
            return [format, files] 
    }
    | set{ input_formats }

    /*
    Files in the fast5 format are converted to pod5 and so are branched out into their respective channels
    */

    input_formats
    | branch { format, files ->
        fast5: format == "fast5"
            return files

        pod5: format == "pod5"
           return files
    }
    | set{ raw_files }

    CONVERT_FAST5_TO_POD5(raw_files.fast5)

    MERGE_POD5(raw_files.pod5)
    | mix(CONVERT_FAST5_TO_POD5.out.pod5_ch) //mix in files if there are only fast5's
    | BASECALL
    | DEMUX
    | flatten
    | map{ long_read_bam -> 
        def meta = [:]
        meta.barcode_kit = params.barcode_kit_name
        meta.barcode = "${ long_read_bam.simpleName.contains("barcode") ? long_read_bam.simpleName.split("barcode")[-1] : long_read_bam.simpleName.split("_")[-1] }" //i.e. when simpleName = unclassified
        tuple(meta, long_read_bam)
    }
    | set{ bam_ch }

    DORADO_SUMMARY(BASECALL.out.called_channel)
    | PYCOQC

    if (params.additional_metadata) {
        bam_ch
        | map { meta, reads -> ["${meta.barcode_kit}_${meta.barcode}", meta, reads]}
        | set { bam_by_barcode_ch }

        Channel.fromPath(params.additional_metadata)
        | ifEmpty {exit 1, "${params.additional_metadata} appears to be an empty file!"}
        | splitCsv(header:true, sep:',')
        | map { meta -> ["${params.barcode_kit_name}_${meta.barcode}", meta] }
        | join(bam_by_barcode_ch)
        | map { barcodekit_barcode, meta1, meta2, reads -> [meta1 + meta2, reads] }
        | set { bam_ch }

    } else {

        bam_ch
        | set { bam_ch }
    }
    
    if (params.read_format == "fastq") {
        CONVERT_TO_FASTQ(bam_ch)
        | set { read_ch }

    } else {
        PUBLISH_BAMS(bam_ch)

        bam_ch
        | set { read_ch }
    }

    emit:
    read_ch
}
