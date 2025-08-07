include { CHECKM2 as PRE_CHECKM2;
          CHECKM2                          } from './modules/checkm2.nf'
include { GTDBTK                           } from './modules/gtdbtk.nf'
include { GUNC as PRE_GUNC;
          GUNC                             } from './modules/gunc.nf'
include { MDMCLEANER                       } from './modules/mdmcleaner.nf'
include { SEQKIT                           } from './modules/seqkit.nf'
include { REPORT                           } from './modules/reporting.nf'
include { FILTER_METADATA as FILTER_REPORT } from '../mixed_input/modules/filter_metadata.nf'
include { FILTER_FASTAS;
          PUBLISH_RESULTS                  } from './modules/filtering.nf'

// Helper functions

def ensureList(Object maybeListLike) {
    if (maybeListLike instanceof Collection) {
        return maybeListLike as List
    } else {
        return [maybeListLike]
    }
}

// Workflows

workflow QC_MAGS {
    take:
    fastas  // [meta, [1.fasta, 2.fasta, ...]]

    main:
    fastas
    | (PRE_CHECKM2 & PRE_GUNC & GTDBTK)

    fastas
    | MDMCLEANER
    | map { meta, fastas ->
        def size = ensureList(fastas).size()
        def group_key = groupKey(meta, size)
        [group_key, meta, fastas]
    }
    | transpose
    | SEQKIT
    | groupTuple
    | map { group_key, meta, fasta_list ->
        [ meta.first(), fasta_list ]
    }
    | set { postqc_fastas }

    postqc_fastas
    | (CHECKM2 & GUNC)

    PRE_CHECKM2.out.results
    | join(PRE_GUNC.out.results)
    | join(CHECKM2.out.results)
    | join(GUNC.out.results)
    | join(GTDBTK.out.results)
    | combine(Channel.fromPath(params.report_config))
    | REPORT

    if (params.filter_manifest) {
        Path filter_manifest = file(params.filter_manifest, checkIfExists: true)
        String fasta_ext = params.fasta_ext.replaceAll(/^\./, '')

        FILTER_REPORT(
            REPORT.out.report,
            filter_manifest,
            "", // select all columns
            false,  // don't remove header
            "" // don't drop duplicates
        )

        FILTER_REPORT.out.filtered_metadata
        | splitCsv(header: true, sep: '\t')
        | map { meta, row -> [meta, "${row.postqc_genome_name}.${fasta_ext}"] }
        | groupTuple
        | set { filtered_fasta_filenames }

        filtered_fasta_filenames
        | join(postqc_fastas)
        | FILTER_FASTAS

        PUBLISH_RESULTS(
            FILTER_REPORT.out.filtered_metadata,
            "${params.outdir}/pass/report"
        )

        // postqc_fastas
        // | join(filtered_fasta_filenames)
        // | map { meta, postqc_fastas_list, filtered_fasta_filenames_list -> 
        //     def postqc_filtered_fastas_list = []
        //     postqc_fastas_list.each { fasta -> 
        //         if (filtered_fasta_filenames_list.contains(fasta.getBaseName())) {
        //             postqc_filtered_fastas_list << fasta
        //         }
        //     }
        //     [meta, postqc_filtered_fastas_list]
        // }
        // | set { filtered_postqc_fastas }
    }
}
