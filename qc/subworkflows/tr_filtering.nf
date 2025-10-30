include { FASTQ2FASTA        } from "../modules/fastq2fasta.nf"
include { TRF                } from "../modules/trf.nf"
include { RMREPEATFROMFASTQ  } from "../modules/rmRepeatFromFq.nf"


workflow TR_FILTERING { 

    /*
    -----------------------------------------------------------------
    Remove tandem repeats from fastq files
    -----------------------------------------------------------------

    */

    take:
    reads_ch // tuple (meta, read_1, read_2)

    main:        // run trf


    // preapare fastq channel to be join by id
    reads_ch.map{meta, fq_1, fq_2 -> 
        tuple (meta.ID, meta, [fq_1, fq_2])
        }
        .set {fqs_ch}
                    
    // convert
    FASTQ2FASTA(reads_ch)
    FASTQ2FASTA.out // tuple (meta, fasta_1, fasta_2)
        | set {trf_in_ch}
    
    TRF(trf_in_ch)
    TRF.out.fasta_trfs // tuple (meta, trf_out_1, trf_out_2)
        | map {meta, trf_out_1, trf_out_2 -> 
            tuple (meta.ID,[trf_out_1, trf_out_2])
            }
        | set {trf_ch} // tuple (meta.id, [trfs_out])

    fqs_ch // tuple (meta.id, meta, [fqs])
        | join(trf_ch) // tuple (meta.id, meta, [fqs], [trfs_out])
        | map {id, meta, fqs, trfs -> tuple(meta, fqs[0], fqs[1], trfs[0], trfs[1])}
        | set {rmTRFfromFq_In_ch}
    RMREPEATFROMFASTQ(rmTRFfromFq_In_ch)
    RMREPEATFROMFASTQ.out.fastqs
        | set {trf_Out_ch} // tuple (meta, trf_fq_gz_1, trf_fq_gz_2)

    emit:
        trf_Out_ch

}
