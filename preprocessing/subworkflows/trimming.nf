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

    if (params.keep_unpaired) {
        unpaired = TRIMMOMATIC.out.trimmed_unpaired_fastqs
    } else {
        unpaired = Channel.empty()
    }
    
    emit:
    paired_ch   = TRIMMOMATIC.out.trimmed_fastqs // tuple (meta, reads_trimmed_1.fq, reads_trimmed_2.fq)
    stats_ch    = COLLATE_STATS_TRIMMOMATIC.out.trimming_stats_ch
    unpaired_ch = unpaired

}