include {TRIMMOMATIC} from "../modules/trimmomatic.nf"
include {FASTQ2FASTA} from "../modules/fastq2fasta.nf"
include {TRF} from "../modules/trf.nf"
include {RMREPEATFROMFASTQ} from "../modules/rmRepeatFromFq.nf"
include {SRA_HUMAN_SCRUBBER} from "../modules/scrubber.nf"
include {COMPRESS_READS; RENAME_READS} from "../modules/helper_processes.nf"
include {SEQTK_MERGEPE; SEQTK_SPLIT} from "../modules/seqtk.nf"

workflow PREPROCESSING {
    /*
    -----------------------------------------------------------------
    Preprocessing fastq files
    -----------------------------------------------------------------

    */

    take:

        reads_ch // tuple (meta, read_1, read_2)

    main:
        // run trimmomatic
        if (params.run_trimmomatic){
            TRIMMOMATIC(reads_ch)
            trimmomatic_Out_ch = TRIMMOMATIC.out.paired_channel // tuple (meta, read_trim_1, read_trim_2)
        }

        // run trf
        if (params.run_trf){
            if (params.run_trimmomatic){
                fastq2fasta_in_ch = trimmomatic_Out_ch
            } else {
                fastq2fasta_in_ch = reads_ch
            }
            
            // preapare fastq channel to be join by id
            fastq2fasta_in_ch.map{meta, fq_1, fq_2 -> 
                tuple (meta.id, meta, [fq_1, fq_2])
                }
                .set {fqs_ch}
            
            // convert
            FASTQ2FASTA(fastq2fasta_in_ch)
            FASTQ2FASTA.out // tuple (meta, fasta_1, fasta_2)
                | set {trf_in_ch}
            
            TRF(trf_in_ch)
            TRF.out.paired_trf // tuple (meta, trf_out_1, trf_out_2)
                | map {meta, trf_out_1, trf_out_2 -> 
                    tuple (meta.id,[trf_out_1, trf_out_2])
                 }
                | set {trf_ch} // tuple (meta.id, [trfs_out])

            fqs_ch // tuple (meta.id, meta, [fqs])
                | join(trf_ch) // tuple (meta.id, meta, [fqs], [trfs_out])
                | map {id, meta, fqs, trfs -> tuple(meta, fqs[0], fqs[1], trfs[0], trfs[1])}
                | set {rmTRFfromFq_In_ch}
            RMREPEATFROMFASTQ(rmTRFfromFq_In_ch)
            RMREPEATFROMFASTQ.out.fastqs
                | set {trf_Out_ch} // tuple (meta, trf_fq_1, trf_fq_2)
        }

        // run human-sra-scrubber
        if (params.run_hrr){

            // if only scrubber is on
            if ((!params.run_trimmomatic) && (!params.run_trf)){
                hrr_In_ch = reads_ch
            }

            // if trimmomatic and no trf
            if ((params.run_trimmomatic) && (!params.run_trf)){
                hrr_In_ch = trimmomatic_Out_ch
            }

            // if trf is true
            if (params.run_trf){
                hrr_In_ch = trf_Out_ch
            }
            // generate interleaved fq for scrubber
            SEQTK_MERGEPE(hrr_In_ch)

            // run hrr
            SRA_HUMAN_SCRUBBER(SEQTK_MERGEPE.out)

            // deinterleaved fastq files
            SEQTK_SPLIT(SRA_HUMAN_SCRUBBER.out, "clean")

            scrubber_Out_ch = SEQTK_SPLIT.out // tuple(meta, reads_clean_1, reads_clean_2)

        }

        // setup output channel
        // if trimmomatics on, trf off and scrubber off
        if ((params.run_trimmomatic) && (!params.run_trf) && (!params.run_hrr)){
            out_ch = trimmomatic_Out_ch // tuple (meta, reads_trim_1, reads_trim_2)
        }
        // if trf on, scrubber off
        if ((params.run_trf) && (!params.run_hrr)){
            out_ch = trf_Out_ch // tuple (meta, trf_fq_1, trf_fq_2)
        }
        // if scrubber on
        if (params.run_hrr){
            out_ch = scrubber_Out_ch // tuple (meta, reads_clean_1, reads_clean_2)
        }

        // publish compressed clean reads
        if (params.publish_clean_reads){
            if (params.compress_clean_reads){
                COMPRESS_READS(out_ch)
            } else {
                RENAME_READS(out_ch)
            }
        }
        // remove unpaired sequences at the end of the process

    emit:
        out_ch // tuple (meta, reads_1, reads_2)
}

def parse_mnf_meta(preprocessing_mnf) {
    // consensus_mnf <Channel.fromPath()>
    println(preprocessing_mnf)
    def mnf_ch =  Channel.fromPath(preprocessing_mnf)
                        | splitCsv(header: true, sep: ',')
                        | map {row -> 
                            // set meta
                            meta = [
                                id: row.sample_id,
                                sample_id: row.sample_id
                                ]

                            // set files
                            reads = [row.reads_1, row.reads_2]

                            // declare channel shape
                            tuple(meta, reads)
                        }
    return mnf_ch // tuple(index, [fastq_pairs])
}

// run preprocessing in isolation
workflow {
    mnf_ch = parse_mnf_meta(params.preprocessing_mnf)
    PREPROCESSING(mnf_ch)
}

def check_preprocessing_params(){
    /*
    -----------------------------------------------------------------
    Checks for necessary parameters and validates paths to ensure 
    they exist. Logs errors if any required parameters are missing.
    -----------------------------------------------------------------

    - **Output**: Number of errors encountered during the checks.

    -----------------------------------------------------------------

    */
    def errors = 0
    // was the kraken database provided?
    if (params.adapter_fasta == null){
        log.error("No adapter_fasta path provided")
        errors +=1
    }

    // if yes, is it a file which exists?
    if (params.adapter_fasta){
        adapter_fasta = file(params.adapter_fasta)
        if (!adapter_fasta.exists()){
            log.error("The adapter_fasta provided (${params.adapter_fasta}) does not exist.")
            errors += 1
        }
    }

    // check switchs

    if ((params.run_trimmomatic == false) && (params.run_trf == false) && (params.run_hrr == false)){
        log.error("All PREPROCESSING process switchs are off (run_trimmomatic = ${params.run_trimmomatic}; run_trf = ${params.run_trf} ; run_hrr = ${params.run_hrr}).")
    }

    return errors
}