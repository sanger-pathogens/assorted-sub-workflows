include { CHECKM2 as PRE_CHECKM2;
          CHECKM2                          } from './modules/checkm2.nf'
include { GTDBTK                           } from './modules/gtdbtk.nf'
include { GUNC as PRE_GUNC;
          GUNC                             } from './modules/gunc.nf'
include { getQuastThresholds;
          QUAST as PRE_QUAST;
          QUAST                            } from './modules/quast.nf'
include { INPUT_FASTA_CHECK                } from './modules/input_check.nf'
include { CLEAN_CONTIGS                    } from './subworkflows/clean_contigs.nf'
include { REPORT;
          COMBINE_REPORTS;
          SUMMARISE_CONTIG_FILTERING       } from './modules/reporting.nf'
include { FILTER_METADATA as FILTER_REPORT } from '../mixed_input/modules/filter_metadata.nf'
include { FILTER_FASTAS;
          PUBLISH_RESULTS                  } from './modules/filtering.nf'
include { ensureList;
          fastaExt;
          genomeName;
          batchFastas                      } from './modules/helper_functions.nf'

// Workflows

workflow QC_MAGS {
    take:
    fasta_lists  // [meta, [1.fasta, 2.fasta, ...]]

    main:

    // genome -> manifest sample ID, recorded BEFORE batching (a batch ID says
    // nothing about which sample a genome came from)
    fasta_lists
    | flatMap { meta, fasta_list ->
        ensureList(fasta_list).collect { f -> "${meta.ID}\t${genomeName(f)}" }
    }
    | collectFile(name: 'genome_sample_map.tsv', newLine: true, sort: true)
    | set { genome_sample_map }

    if (params.batch_size > 0) {
        log.info "Batch size parameter is set to ${params.batch_size}, will process in batches of this size."
        fastas = batchFastas(fasta_lists)
    } else {
        fastas = fasta_lists
    }

    // Clean and check every input genome once per batch before any QC step reads it (modules/input_check.nf).
    INPUT_FASTA_CHECK(fastas)

    INPUT_FASTA_CHECK.out.fixes
    | collectFile(name: 'input_fasta_fixes.tsv', storeDir: "${params.outdir}/report", sort: true)
    | subscribe { f ->
        def n = f.readLines().count { line -> line.trim() }
        if (n > 0) {
            log.warn "[input check] ${n} fix(es)/warning(s) on input genomes (e.g. GenBank '//' lines removed, " +
                     "non-nucleotide characters); see ${params.outdir}/report/input_fasta_fixes.tsv"
        }
    }

    checked_fastas = INPUT_FASTA_CHECK.out.fastas

    pre_qc = Channel.value("pre_qc")
    
    PRE_CHECKM2(checked_fastas, pre_qc)
    PRE_GUNC(checked_fastas, pre_qc)
    PRE_QUAST(checked_fastas, pre_qc, getQuastThresholds())
    
    CLEAN_CONTIGS(checked_fastas)
    postqc_fastas = CLEAN_CONTIGS.out.fastas

    post_qc = Channel.value("post_qc")

    QUAST(postqc_fastas, post_qc, getQuastThresholds())
    CHECKM2(postqc_fastas, post_qc)
    GUNC(postqc_fastas, post_qc)

    GTDBTK(postqc_fastas, post_qc)

    PRE_QUAST.out.results
        | join(PRE_CHECKM2.out.results)
        | join(PRE_GUNC.out.results)
        | join(QUAST.out.results)
        | join(CHECKM2.out.results)
        | join(GUNC.out.results)
        | join(GTDBTK.out.results)
        | combine(Channel.fromPath(params.report_config))
        | REPORT

    COMBINE_REPORTS(
        REPORT.out.report.map { meta, report -> report }.collect(),
        genome_sample_map
    )
        
    PRE_QUAST.out.results
        | join(QUAST.out.results)
        | SUMMARISE_CONTIG_FILTERING

    // A batch is finished once its report and contig-filtering summary exist (and, with autoqc, its pass outputs)
    batch_finished = REPORT.out.report
        .join(SUMMARISE_CONTIG_FILTERING.out.summary)
        .map { meta, report, summary -> [meta, 'done'] }

    if (params.autoqc_config) {
        // Check if user-supplied config or opted for default config
        Path autoqc_config = params.autoqc_config == "default"
            ? file("${projectDir}/assorted-sub-workflows/qc_mags/assets/autoqc_config.tsv", checkIfExists: true)
            : file(params.autoqc_config, checkIfExists: true)

        String fasta_ext = fastaExt()

        FILTER_REPORT(
            REPORT.out.report,
            autoqc_config,
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

        batch_finished = batch_finished
            .join(PUBLISH_RESULTS.out)
            .join(FILTER_FASTAS.out.filtered_fastas, remainder: true)
            .map { meta, _finished, _published, _passed -> [meta, 'done'] }
    }

    if (!params.skip_cleanup) {
        INPUT_FASTA_CHECK.out.fastas
        | mix(
            PRE_CHECKM2.out.results,
            PRE_GUNC.out.results,
            PRE_QUAST.out.results,
            CLEAN_CONTIGS.out.mdmcleaner,
            postqc_fastas,
            QUAST.out.results,
            CHECKM2.out.results,
            GUNC.out.results,
            GTDBTK.out.results
        )
        | combine(batch_finished, by: 0)
        | subscribe { _meta, outputs, _finished -> NextflowTool.cleanWork(outputs, workflow.workDir, log) }
    }
}
