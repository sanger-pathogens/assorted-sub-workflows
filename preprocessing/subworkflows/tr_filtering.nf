include { FASTQ2FASTA              } from "../modules/fastq2fasta.nf"
include { TRF                      } from "../modules/trf.nf"
include { RM_REPEAT_FROM_PAIRED;
          RM_REPEAT_FROM_UNPAIRED  } from "../modules/rmRepeatFromFq.nf"


workflow TR_FILTERING { 

    /*
    -----------------------------------------------------------------
    Remove tandem repeats from fastq files
    -----------------------------------------------------------------
    */

    take:
    reads_ch // tuple (meta, read_1, read_2)
    unpaired_ch // tuple (meta, unpaired_read)

    main:

    /* 
        prepare fastq channel to be join by id
    */
    reads_ch
    | flatMap { meta, fq_1, fq_2 ->
        [
            tuple(meta.ID, "1", meta, fq_1),
            tuple(meta.ID, "2", meta, fq_2)
        ]
    }
    | set { flat_paired_ch }

    unpaired_ch
    | map { meta, fq -> tuple(meta.ID, "unpaired", meta, fq) }
    | set { flat_unpaired_ch }

    flat_paired_ch
    | mix(flat_unpaired_ch)
    | set { flat_reads_ch }
                    
    /*
        actual pipeline
    */

    flat_reads_ch
    | FASTQ2FASTA
    | TRF // tuple (meta.ID, read_type, meta, fastq, trf_file)
    | branch { id, read_type, meta, fastq, trf ->
        paired:   read_type == "1" || read_type == "2"
        unpaired: read_type == "unpaired"
    }
    | set { trf_branched }

    trf_branched.paired
    | groupTuple(by: 0, size: 2)
    | map { id, read_types, metas, fastqs, trfs ->
        assert read_types.size() == 2 :
            "Sample ${id} expected R1 and R2 but got: ${read_types}"

        def meta  = metas[0]
        
        def index1 = read_types.indexOf("1")
        def index2 = read_types.indexOf("2")

        tuple(meta, fastqs[index1], fastqs[index2], trfs[index1], trfs[index2])
    }
    | RM_REPEAT_FROM_PAIRED


    trf_branched.unpaired
    | map { id, read_type, meta, fastq, trf -> 
        tuple(meta, fastq, trf) 
    }
    | RM_REPEAT_FROM_UNPAIRED


    emit:
    trf_paired_ch   = RM_REPEAT_FROM_PAIRED.out.fastqs // tuple (meta, read_1, read_2)
    trf_unpaired_ch = RM_REPEAT_FROM_UNPAIRED.out.fastqs // tuple (meta, unpaired_read)
}