include { FASTQC } from '../modules/fastqc.nf'
include { MULTIQC } from '../../reporting/modules/multiqc.nf'

workflow FASTQC_MULTIQC {

    take:
    fastq_path_ch

    main:
    if (!params.skip_fastqc) {
        FASTQC(fastq_path_ch)
        def post_qc_report = false

        FASTQC.out.zip
        | map { meta, zip1, zip2 -> [zip1, zip2] }
        | collect()
        | flatten()
        | set { fastqc_zips }

        MULTIQC(fastqc_zips, post_qc_report)
    }
    
    emit:
    fastqc_report = MULTIQC.out.report
}