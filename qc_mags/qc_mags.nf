include { CHECKM2 as PRE_CHECKM2;
          CHECKM2                   } from '../modules/qc_mags/checkm2.nf'
include { GUNC as PRE_GUNC;
          GUNC                      } from '../modules/qc_mags/gunc.nf'
include { MDMCLEANER                } from '../modules/qc_mags/mdmcleaner.nf'
include { SEQKIT                    } from '../modules/qc_mags/seqkit.nf'
include { REPORT                    } from '../modules/qc_mags/reporting.nf'


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