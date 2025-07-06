include { PUBLISH_RESULTS; SUMMARY } from '../modules/reporting/report.nf'

workflow REPORT_BINS {
    
    take:
    bin_ch

    main:

    PUBLISH_RESULTS(bin_ch)

    bin_ch
    | map { meta, fastas, checkmReport -> [checkmReport]}
    | collect
    | SUMMARY
}