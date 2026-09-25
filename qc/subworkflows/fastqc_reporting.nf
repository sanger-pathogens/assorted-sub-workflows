include { FASTQC } from '../modules/fastqc.nf'
include { MULTIQC } from '../../reporting/modules/multiqc.nf'

workflow FASTQC_MULTIQC {

    take:
    fastq_path_ch

    main:

    if (!params.skip_fastqc) {
        FASTQC(fastq_path_ch)

        // use method-call chaining: in a pipe, `collect()` with parentheses resolves
        // to Groovy's Object.collect() rather than the Nextflow operator
        fastqc_zips = FASTQC.out.zip
            .map { meta, zip1, zip2 -> [zip1, zip2] }
            .flatten()
            .collect()

        MULTIQC(fastqc_zips, []) // no custom multiqc config
        fastqc_report = MULTIQC.out.report
    } else {
        fastqc_report = Channel.value("FastQC skipped")
    }

    emit:
    fastqc_report
}
