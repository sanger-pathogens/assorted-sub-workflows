include { CHECKM2 as PRE_CHECKM2;
          CHECKM2                   } from './modules/checkm2.nf'
include { GUNC as PRE_GUNC;
          GUNC                      } from './modules/gunc.nf'
include { MDMCLEANER                } from './modules/mdmcleaner.nf'
include { SEQKIT                    } from './modules/seqkit.nf'
include { REPORT                    } from './modules/reporting.nf'


workflow QC_MAGS {
    take:
    fasta_directory

    main:
    
    fasta_directory
    | (PRE_CHECKM2 & PRE_GUNC)

    fasta_directory
    | MDMCLEANER
    | SEQKIT
    | (CHECKM2 & GUNC)

    PRE_CHECKM2.out.results
    | join(PRE_GUNC.out.results)
    | join(CHECKM2.out.results)
    | join(GUNC.out.results)
    | REPORT
}