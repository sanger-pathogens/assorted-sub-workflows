include { BMTAGGER                    } from "../modules/bmtagger.nf"
include { FILTER_HOST_READS; 
          GET_HOST_READS              } from '../modules/filter_reads.nf'
include { GENERATE_STATS              } from '../modules/generate_stats.nf'
include { COLLATE_STATS_BMTAGGER      } from '../modules/collate_stats.nf'



workflow READ_REMOVAL {
    /*
    -----------------------------------------------------------------
    Remove reads mapping to unwanted references
    -----------------------------------------------------------------

    */

    take:

        reads_ch // tuple (meta, read_1, read_2)

    main:

        // run hrr
    BMTAGGER(reads_ch)
    | set {bmtagger_out_ch}

    FILTER_HOST_READS(BMTAGGER.out.data_ch, BMTAGGER.out.bmtagger_list_ch)
    | set { filtered_reads_ch }

    GET_HOST_READS(BMTAGGER.out.data_ch, BMTAGGER.out.bmtagger_list_ch)
    | set { host_reads_ch }

    filtered_reads_ch.data_ch
    | join(filtered_reads_ch.cleaned_ch)
    | join(host_reads_ch.host_ch)
    | join(reads_ch)
    | GENERATE_STATS

    COLLATE_STATS_BMTAGGER(GENERATE_STATS.out.stats_ch.collect())

    read_removal_Out_ch = filtered_reads_ch.cleaned_ch // tuple (meta, reads_clean_1.gz, reads_clean_2.gz)   
    
    emit: 
    read_removal_Out_ch
}