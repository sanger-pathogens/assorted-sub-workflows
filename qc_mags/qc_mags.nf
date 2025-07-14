include { CHECKM2 as PRE_CHECKM2;
          CHECKM2                   } from './modules/checkm2.nf'
include { GTDBTK                    } from './modules/gtdbtk.nf'
include { GUNC as PRE_GUNC;
          GUNC                      } from './modules/gunc.nf'
include { MDMCLEANER                } from './modules/mdmcleaner.nf'
include { SEQKIT                    } from './modules/seqkit.nf'
include { BUNDLE_FASTAS             } from './modules/helper_scripts.nf'
include { REPORT                    } from './modules/reporting.nf'


workflow QC_MAGS {
    take:
    fasta_directory

    main:
    fasta_directory
    | (PRE_CHECKM2 & PRE_GUNC & GTDBTK)

    fasta_directory
    | MDMCLEANER
    | map { meta, fasta_list ->
        def size = fasta_list.size()
        def group_key = groupKey(meta, size)
        [group_key, meta, fasta_list]
    }
    | transpose
    | SEQKIT
    | groupTuple
    | map { group_key, meta, fasta_list ->
        [ meta.first(), fasta_list ]
    }
    | BUNDLE_FASTAS
    | (CHECKM2 & GUNC)

    PRE_CHECKM2.out.results
    | join(GTDBTK.out.results)
    | join(PRE_GUNC.out.results)
    | join(CHECKM2.out.results)
    | join(GUNC.out.results)
    | REPORT
}