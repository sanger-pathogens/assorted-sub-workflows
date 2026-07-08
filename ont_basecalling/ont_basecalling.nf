import java.nio.file.Files

include { CONVERT_FAST5_TO_POD5; MERGE_POD5; CHECK_POD5_CHEMISTRY } from './modules/pod5.nf'
include { BASECALL; DEMUX; DORADO_SUMMARY; BASECALL_LEGACY        } from './modules/dorado.nf'
include { PYCOQC                                                  } from './modules/pycoqc.nf'
include { CONVERT_TO_FASTQ; PUBLISH_BAMS                          } from './modules/samtools.nf'

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
    /*
    Files in the fast5 format are converted to pod5 and so are branched out into their respective channels
    */

    raw_read_signal_files
    | branch { meta, files ->
        fast5: meta.format == "fast5"
            return tuple(meta, files)

        pod5: meta.format == "pod5"
            return tuple(meta, files)
    }
    | set{ raw_files }

    CONVERT_FAST5_TO_POD5(raw_files.fast5)

    MERGE_POD5(raw_files.pod5)
    | mix(CONVERT_FAST5_TO_POD5.out.pod5_ch) //mix in files if there are only fast5's
    | CHECK_POD5_CHEMISTRY
    | map { meta, report_json, files ->
        def info = new groovy.json.JsonSlurper().parse(report_json)
        def merged_meta = meta + info
        tuple(merged_meta, files)
    }
    | branch { meta, _pod5 ->
        legacy_0_9_6:  meta.dorado_version == '0.9.6'
        current: true
    }
    | set { ready_for_basecalling_ch }
    
    BASECALL(ready_for_basecalling_ch.current)
    BASECALL_LEGACY(ready_for_basecalling_ch.legacy_0_9_6)

    BASECALL.out.called_channel
    | mix(BASECALL_LEGACY.out.called_channel)
    | set { called_ch }

    if (params.barcode_kit_name) {
        validateBarcodeParams()

        called_ch
        | DEMUX
        | flatten
        | map{ long_read_bam -> 
            def meta = [:]
            meta.barcode_kit = params.barcode_kit_name
            meta.barcode = long_read_bam.simpleName.contains("barcode") ? long_read_bam.simpleName.split("barcode")[-1] : long_read_bam.simpleName.split("_")[-1] //i.e. when simpleName = unclassified
            tuple(meta, long_read_bam)
        }
        | set{ bam_ch }
    } else {
        called_ch
        | set{ bam_ch }
    }

    DORADO_SUMMARY(called_ch)
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
        | map { _barcodekit_barcode, meta1, meta2, reads -> [meta1 + meta2, reads] }
        | set { bam_ch }

    } else {
        bam_ch
        | map { meta, called_bam -> 
            //as we are not performing any meta additon we just generate a name from the input directory
            def input_directory_name = file(params.raw_read_dir).simpleName
            meta.barcode_kit = "not-reclassified"
            meta.barcode = input_directory_name.contains("barcode") ? input_directory_name.split("barcode")[-1] : input_directory_name
            [meta, called_bam]
        }
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

    //cleanup pod5
    BASECALL.out.basecalling_complete_channel
    | mix(BASECALL_LEGACY.out.basecalling_complete_channel)
    | flatten
    | filter(Path)
    | map { pod5File -> NextflowTool.safeDelete(pod5File, workflow.workDir, log) }


    emit:
    read_ch
}
