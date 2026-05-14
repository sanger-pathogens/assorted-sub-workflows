include { BMTAGGER_PAIRED;
          BMTAGGER_SINGLE            } from "../modules/bmtagger.nf"
include { FILTER_HOST_READS;
          FILTER_HOST_READS_SINGLE;
          GET_HOST_READS;
          GET_HOST_READS_SINGLE      } from '../modules/filter_reads.nf'
include { GENERATE_STATS             } from '../modules/generate_stats.nf'
include { COLLATE_STATS_BMTAGGER     } from '../modules/collate_stats.nf'



workflow HOST_READ_REMOVAL {
    /*
    -----------------------------------------------------------------
    Remove reads mapping to unwanted references
    -----------------------------------------------------------------

    */

    take:

    reads_ch // tuple (meta, read_1, read_2)
    unpaired_reads_ch // tuple (meta, unpaired_read)

    main:

    // run hrr
    BMTAGGER_PAIRED(reads_ch)
    | GET_HOST_READS & FILTER_HOST_READS

    BMTAGGER_SINGLE(unpaired_reads_ch)
    | GET_HOST_READS_SINGLE & FILTER_HOST_READS_SINGLE

    FILTER_HOST_READS.out.data_ch
    | join(FILTER_HOST_READS.out.cleaned_ch)
    | join(GET_HOST_READS.out.host_ch)
    | join(reads_ch)
    | GENERATE_STATS

    COLLATE_STATS_BMTAGGER(GENERATE_STATS.out.stats_ch.collect())
    
    emit: 
    paired_ch                    = FILTER_HOST_READS.out.cleaned_ch
    unpaired_ch                  = FILTER_HOST_READS_SINGLE.out.cleaned_ch
    collated_host_reads_stats_ch = COLLATE_STATS_BMTAGGER.out.host_reads_stats_ch
}