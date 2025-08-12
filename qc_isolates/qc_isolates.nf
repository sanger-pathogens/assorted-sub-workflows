include { CHECKM2                   } from './modules/checkm2.nf'
include { GTDBTK                    } from './modules/gtdbtk.nf'
include { GUNC                      } from './modules/gunc.nf'
include { REPORT                    } from './modules/reporting.nf'

workflow QC_ISOLATES {
    take:
    fastas

    main:
    fastas
    | (CHECKM2 & GUNC & GTDBTK)
 
    CHECKM2.out.results
    | join(GUNC.out.results)
    | join(GTDBTK.out.results)
    | combine(Channel.fromPath(params.report_config))
    | REPORT
}