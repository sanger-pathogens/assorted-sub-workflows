include { CHECKM2                   } from './modules/checkm2.nf'
include { GTDBTK                    } from './modules/gtdbtk.nf'
include { GUNC                      } from './modules/gunc.nf'
include { QUAST;                         
          QUAST_SUMMARY             } from './modules/quast.nf'
include { REPORT                    } from './modules/reporting.nf'

workflow QC_ISOLATES {
    take:
    fastas

    main:
    fastas
    | (GTDBTK & QUAST)

    QUAST.out.results | QUAST_SUMMARY

    fastas
    | (SEQKIT)

    SEQKIT.out.results
    | (CHECKM2 & GUNC)
 
    CHECKM2.out.results
    | join(GUNC.out.results)
    | join(GTDBTK.out.results)
    | join(QUAST_SUMMARY.out.results)
    | combine(Channel.fromPath(params.report_config))
    | REPORT
}