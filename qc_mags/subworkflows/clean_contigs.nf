include { MDMCLEANER                      } from '../modules/mdmcleaner.nf'
include { SEQKIT                          } from '../modules/seqkit.nf'
include { ensureList;
          genomeName;
          genomeType;
          genomeTypeMap                   } from '../modules/helper_functions.nf'

workflow CLEAN_CONTIGS {
    take:
    fastas  // [meta, [1.fasta, 2.fasta, ...]]

    main:
    // Split each batch into one genome per task, keyed so the batch can be regrouped once all its genomes are done
    fastas
    | map { meta, fasta_list ->
        def fasta_list_ = ensureList(fasta_list)
        [groupKey(meta, fasta_list_.size()), meta, fasta_list_]
    }
    | transpose
    | set { genomes }

    // MDMCLEANER removes contaminant contigs from MAGs (binning artefacts). Isolates are pure cultures:
    // cleaning them would only risk stripping real accessory content, so they go straight to SEQKIT, which
    // gives both genome types the same cleaned_<name>_filtered_kept_contigs file name.
    def gtypes = genomeTypeMap()
    genomes
    | branch { group_key, meta, fasta ->
        isolate: genomeType(gtypes, genomeName(fasta)) == 'isolate'
        mag: true
    }
    | set { by_type }

    MDMCLEANER(by_type.mag)

    MDMCLEANER.out.results
    | mix(by_type.isolate)
    | SEQKIT
    | groupTuple(remainder: true)
    | map { group_key, meta, fasta_list ->
        [ meta.first(), fasta_list ]
    }
    | set { cleaned_fastas }

    // Audit: MAGs sent to MDMCLEANER with no cleaned output failed and were dropped by errorStrategy 'ignore'
    by_type.mag
    | map { group_key, meta, fasta -> [[meta.ID, genomeName(fasta)], fasta] }
    | join(
        MDMCLEANER.out.results
        | map { group_key, meta, fasta -> [[meta.ID, genomeName(fasta).replaceAll(/_filtered_kept_contigs$/, '')], fasta] },
        remainder: true
    )
    | filter { key, sent, cleaned -> cleaned == null }
    | map { key, sent, cleaned -> "${key.join('\t')}\n" }
    | concat(Channel.of(''))  // so the file is written even when nothing failed
    | collectFile(name: 'mdmcleaner_skipped.tsv', seed: "batch\tgenome\n", storeDir: "${params.outdir}/report")
    | subscribe { f ->
        def n = f.countLines() - 1
        if (n > 0) {
            log.warn "[MDMCLEANER audit] ${n} genome(s) failed MDMCLEANER and have NO post-QC result " +
                     "(dropped by errorStrategy 'ignore'); see ${params.outdir}/report/mdmcleaner_skipped.tsv"
        }
    }

    emit:
    fastas = cleaned_fastas  // [meta, [cleaned_1.fasta, cleaned_2.fasta, ...]]
    mdmcleaner = MDMCLEANER.out.results.map { group_key, meta, fasta -> [meta, fasta] }  // per genome, for cleanup
}
