include { TRIMMOMATIC                    } from "../modules/trimmomatic.nf"
include { COLLATE_STATS_TRIMMOMATIC      } from '../modules/collate_stats.nf'



workflow TRIMMING {
    /*
    -----------------------------------------------------------------
    Remove reads mapping to unwanted references
    -----------------------------------------------------------------

    */

    take:

        reads_ch // tuple (meta, read_1, read_2)

    main:

        TRIMMOMATIC(reads_ch)   

        COLLATE_STATS_TRIMMOMATIC(TRIMMOMATIC.out.trimmomatic_stats.collect())
    
    emit: 

        TRIMMOMATIC.out.trimmed_fastqs // tuple (meta, reads_trimmed_1.fq, reads_trimmed_2.fq)
}