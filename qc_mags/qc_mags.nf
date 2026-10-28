include { CHECKM2 as PRE_CHECKM2;
          CHECKM2                          } from './modules/checkm2.nf'
include { GTDBTK                           } from './modules/gtdbtk.nf'
include { GUNC as PRE_GUNC;
          GUNC                             } from './modules/gunc.nf'
include { getQuastThresholds;
          QUAST as PRE_QUAST;
          QUAST                            } from './modules/quast.nf'
include { INPUT_FASTA_CHECK                } from './modules/input_check.nf'
include { MDMCLEANER;
          MDMCLEANER_SKIP;
          MDMCLEANER_AUDIT                 } from './modules/mdmcleaner.nf'
include { SEQKIT                           } from './modules/seqkit.nf'
include { REPORT;
          COMBINE_REPORTS;
          SUMMARISE_CONTIG_FILTERING       } from './modules/reporting.nf'
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

// Genome type per genome (preqc name -> 'mag'|'isolate') for MIXED runs, from --genome_types (curated
// provenance); genomes not listed use --genome_type. Anything that is not exactly 'isolate' is a MAG.
def genomeTypeMap() {
    if (!params.genome_types) { return [:] }
    def rows = file(params.genome_types, checkIfExists: true).readLines()
    def header = rows[0].split('\t').collect { it.trim() }
    def gi = header.indexOf('genome')
    def ti = header.indexOf('genome_type')
    if (gi < 0 || ti < 0) { error "--genome_types ${params.genome_types}: needs columns 'genome' and 'genome_type'" }
    rows.drop(1).findAll { it.trim() }.collectEntries { l ->
        def c = l.split('\t', -1)
        [(c[gi].trim()): c[ti].trim().toLowerCase()]
    }
}

// Workflows

workflow QC_MAGS {
    take:
    fasta_lists  // [meta, [1.fasta, 2.fasta, ...]]

    main:

    // genome -> manifest sample ID, recorded BEFORE batching (a batch ID says
    // nothing about which sample a genome came from). Keyed on the fasta name
    // without the extension, which is what the reports call preqc_genome_name.
    String sample_map_ext = ".${params.fasta_ext.replaceAll(/^\./, '')}"
    fasta_lists
    | flatMap { meta, fasta_list ->
        ensureList(fasta_list).collect { f ->
            def name = f.name.endsWith(sample_map_ext) ? f.name[0..<(f.name.length() - sample_map_ext.length())] : f.baseName
            "${meta.ID}\t${name}"
        }
    }
    | collectFile(name: 'genome_sample_map.tsv', newLine: true, sort: true)
    | set { genome_sample_map }

    if (params.batch_size > 0) {
        log.info "Batch size parameter is set to ${params.batch_size}, will process in batches of this size."

        def batchIndex = 0

        fasta_lists
        | map { meta, fasta_list ->
            fasta_list
        }
        | flatten()
        | buffer( size: params.batch_size, remainder: true )
        | map { batch -> 
            def meta = [:]
            meta.ID = "batch_${++batchIndex}"
            [ meta, batch ]
        }
        | set { fastas }
    } else {
        fasta_lists
        | set { fastas }
    }

    // Clean and check every input genome once per batch before any QC step reads it (modules/input_check.nf).
    INPUT_FASTA_CHECK(fastas)
    
    INPUT_FASTA_CHECK.out.fixes
    | collectFile(name: 'input_fasta_fixes.tsv', storeDir: "${params.outdir}/report", sort: true)
    | map { f -> f.readLines().findAll { line -> line.trim() }.size() }
    | subscribe { n ->
        if (n > 0) {
            log.warn "[input check] ${n} fix(es)/warning(s) on input genomes (e.g. GenBank '//' lines removed, " +
                     "non-nucleotide characters); see ${params.outdir}/report/input_fasta_fixes.tsv"
        }
    }

    INPUT_FASTA_CHECK.out.fastas 
    | set { checked_fastas }

    pre_qc = Channel.value("pre_qc")
    
    PRE_CHECKM2(checked_fastas, pre_qc)
    PRE_GUNC(checked_fastas, pre_qc)
    PRE_QUAST(checked_fastas, pre_qc, getQuastThresholds())
    
    checked_fastas
    | map { meta, fastas ->
        def size = ensureList(fastas).size()
        def group_key = groupKey(meta, size)
        [group_key, meta, fastas]
    }
    | transpose
    | set { mdmcleaner_input }

    // MDMCLEANER removes contaminant contigs from MAGs (binning artefacts). Isolates are pure cultures:
    // cleaning them would only risk stripping real accessory content, so they skip it. MDMCLEANER_SKIP gives
    // them the same <name>_filtered_kept_contigs file name, so SEQKIT, post-QC tools and the report joins
    // are identical for both genome types.
    def gtypes = genomeTypeMap()
    String run_type = (params.genome_type == 'isolate') ? 'isolate' : 'mag'
    String type_ext = ".${params.fasta_ext.replaceAll(/^\./, '')}"
    mdmcleaner_input
    | branch { group_key, meta, fasta ->
        def name = fasta.name.endsWith(type_ext) ? fasta.name[0..<(fasta.name.length() - type_ext.length())] : fasta.baseName
        isolate: (gtypes.containsKey(name) ? gtypes[name] : run_type) == 'isolate'
        mag: true
    }
    | set { by_type }

    MDMCLEANER(by_type.mag)
    MDMCLEANER_SKIP(by_type.isolate)

    MDMCLEANER.out.results
    | mix(MDMCLEANER_SKIP.out.results)
    | SEQKIT
    | groupTuple(remainder: true)
    | map { group_key, meta, fasta_list ->
        [ meta.first(), fasta_list ]
    }
    | set { postqc_fastas }

    String mdm_ext = ".${params.fasta_ext.replaceAll(/^\./, '')}"
    by_type.mag   // isolates never go to MDMCLEANER, so they are not "sent"
    | map { group_key, meta, fasta ->
        def name = fasta.name.endsWith(mdm_ext) ? fasta.name[0..<(fasta.name.length() - mdm_ext.length())] : fasta.baseName
        "${meta.ID}\t${name}"
    }
    | concat(Channel.of(''))
    | collectFile(name: 'mdmcleaner_sent.tsv', newLine: true)
    | set { mdmcleaner_sent }

    MDMCLEANER.out.results
    | map { group_key, meta, fasta -> "${meta.ID}\t${fasta.name.replaceAll(/_filtered_kept_contigs\..*$/, '')}" }
    | concat(Channel.of(''))
    | collectFile(name: 'mdmcleaner_cleaned.tsv', newLine: true)
    | set { mdmcleaner_cleaned }

    MDMCLEANER_AUDIT(mdmcleaner_sent, mdmcleaner_cleaned)
    MDMCLEANER_AUDIT.out.n_skipped
    | subscribe { n ->
        if (n.trim().isInteger() && n.trim().toInteger() > 0) {
            log.warn "[MDMCLEANER audit] ${n.trim()} genome(s) failed MDMCLEANER and have NO post-QC result " +
                     "(dropped by errorStrategy 'ignore'); see ${params.outdir}/report/mdmcleaner_skipped.tsv"
        }
    }

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

    if (params.autoqc_config) {
        // Check if user-supplied config or opted for default config
        Path autoqc_config = params.autoqc_config == "default"
            ? file("${projectDir}/assorted-sub-workflows/qc_mags/assets/autoqc_config.tsv", checkIfExists: true)
            : file(params.autoqc_config, checkIfExists: true)

        String fasta_ext = params.fasta_ext.replaceAll(/^\./, '')

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
    }
}
