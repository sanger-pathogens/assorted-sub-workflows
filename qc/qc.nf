#!/usr/bin/env nextflow
include { BOWTIE2_HOST_REMOVAL } from './modules/host_removal.nf'
include { FASTP
          COLLATE_FASTP_STATS  } from './modules/fastp.nf'
include { SYLPH_SKETCH
          SYLPH_PROFILE
          SYLPHTAX_TAXPROF     } from './modules/sylph.nf'

//switch is banned :(
def needs_host_removal(meta) {
    if (params.host_removal == "all") {
        return true
    }
    if (params.host_removal == "none") {
        return false
    }
    // "auto": reads from iRODS have already been through NPG's human read removal
    return meta.input_source != "IRODS"
}

workflow QC {
    take:
    reads_ch // meta, read_1, read_2

    main:
    // items are [meta, read_1, read_2] (paired) or [meta, reads] (ONT), so only meta is read here
    reads_ch.branch{ reads ->
        // PUBLIC_PAIRED: paired reads from ENA (see MIXED_INPUT), QC'd like Illumina
        illumina: reads[0].Platform in ["ILLUMINA", "PUBLIC_PAIRED"]

        ONT: reads[0].Platform == "ONT"

        other: true
    }
    | set { reads_ch }

    // only Illumina reads are QC'd for now, so pass everything else through
    reads_ch.ONT
    | mix(reads_ch.other)
    | set { not_qcd_ch }

    not_qcd_ch.subscribe { reads ->
        log.info "${reads[0].ID}: passed through without QC, as only Illumina reads are QC'd (Platform: ${reads[0].Platform})"
    }

    QC_ILLUMINA(reads_ch.illumina)

    emit:
    preprocessed_reads_ch = QC_ILLUMINA.out.preprocessed_reads_ch.mix(not_qcd_ch)
    fastp_stats = QC_ILLUMINA.out.fastp_stats
    host_removal_stats = QC_ILLUMINA.out.host_removal_stats
    sylph_report = QC_ILLUMINA.out.sylph_report
    sylphtax_mpa_report = QC_ILLUMINA.out.sylphtax_mpa_report
}

workflow QC_ILLUMINA {
    /*
    -----------------------------------------------------------------
    QC of paired short reads:
    bowtie2 human read removal (see host_removal) -> fastp trimming -> sylph taxonomic profiling
    -----------------------------------------------------------------
    */

    take:
    reads_ch // meta, read_1, read_2

    main:
    if (params.run_qc) {
        if (!(params.host_removal in ["auto", "all", "none"])) {
            error "--host_removal must be one of auto, all or none (got '${params.host_removal}')"
        }

        reads_ch.branch { meta, read_1, read_2 ->
            host_removal: needs_host_removal(meta)
            no_host_removal: true
        }
        | set { reads_by_host_removal_ch }

        BOWTIE2_HOST_REMOVAL(reads_by_host_removal_ch.host_removal)

        BOWTIE2_HOST_REMOVAL.out.stats
        | collectFile(name: "host_removal_statistics.tsv", keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/preprocessing_summary_stats")
        | set { host_removal_stats }

        BOWTIE2_HOST_REMOVAL.out.reads
        | mix(reads_by_host_removal_ch.no_host_removal)
        | FASTP

        FASTP.out.json.collect()
        | COLLATE_FASTP_STATS
        | set { fastp_stats }

        preprocessed_reads_ch = FASTP.out.reads

        if (params.sylph_profile) {
            SYLPH(preprocessed_reads_ch)
            sylph_report = SYLPH.out.sylph_report
            sylphtax_mpa_report = SYLPH.out.sylphtax_mpa_report
        } else {
            sylph_report = Channel.empty()
            sylphtax_mpa_report = Channel.empty()
        }

    } else {
        preprocessed_reads_ch = reads_ch
        fastp_stats = Channel.empty()
        host_removal_stats = Channel.empty()
        sylph_report = Channel.empty()
        sylphtax_mpa_report = Channel.empty()
    }

    emit:
    preprocessed_reads_ch
    fastp_stats
    host_removal_stats
    sylph_report
    sylphtax_mpa_report
}

workflow SYLPH {
    take:
    reads_ch // meta, read_1, read_2

    main:
    SYLPH_SKETCH(reads_ch)

    // batches are sorted by name so the same samples give the same batches and -resume reuses them
    SYLPH_SKETCH.out.sketch
    | map { meta, sketch -> sketch }
    | toSortedList { a, b -> a.name <=> b.name }
    | flatMap { sketches -> sketches.collate(params.sylph_profile_batch_size).withIndex().collect { batch, index -> [ index, batch ] } }
    | SYLPH_PROFILE

    SYLPH_PROFILE.out.batch_report
    | combine(channel.fromPath(params.sylph_tax_metadata, checkIfExists: true))
    | SYLPHTAX_TAXPROF

    // per-sample reports are named after the sample ID, so pair them back up with their meta
    SYLPH_SKETCH.out.sketch
    | map { meta, sketch -> [ meta.ID, meta ] }
    | set { meta_by_id_ch }

    SYLPH_PROFILE.out.sylph_reports
    | flatten
    | map { report -> [ report.name - ~/_sylph_profile\.tsv$/, report ] }
    | join(meta_by_id_ch)
    | map { id, report, meta -> [ meta, report ] }
    | set { sylph_report }

    SYLPHTAX_TAXPROF.out.sylphtax_mpa_reports
    | flatten
    | map { report -> [ report.baseName, report ] }
    | join(meta_by_id_ch)
    | map { id, report, meta -> [ meta, report ] }
    | set { sylphtax_mpa_report }

    emit:
    sylph_report
    sylphtax_mpa_report
}