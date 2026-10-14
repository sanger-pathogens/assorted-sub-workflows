include { BEEFEATER_SEARCH;
          BEEFEATER_GET                 } from '../modules/beefeater.nf'
include { COLLATE_FASTQ                 } from '../modules/samtools.nf'
include { METADATA as METADATA_QUERIED  } from '../modules/metadata_save.nf'
include { PUBLISH_FASTQ
          PUBLISH_UNBASECALLED          } from '../modules/publish.nf'

/*
Fail on records that share an ID. beefeater names each download by ID, so records sharing one
would overwrite each other's file, both be paired with whichever copy landed, and collide on publish.
*/
def check_unique_ids(records) {
    def duplicates = records.groupBy { it.ID }.findAll { id, group -> group.size() > 1 }
    if (duplicates) {
        def details = duplicates.collect { id, group -> "  ${id}: ${group*.irods_path.join(', ')}" }.join('\n')
        error "${duplicates.size()} ID(s) are shared by more than one search record, so their downloads would overwrite each other:\n${details}"
    }
}

/*
Group search records into json batch files of batch_size, named <prefix>_batch_<n>.json.
Records are sorted by ID first, as beefeater emits them in download-completion order,
so the same records give the same batches and -resume reuses completed downloads.
*/
def batch_records(records_ch, batch_size, prefix) {
    return records_ch
        .toSortedList { a, b -> a.ID <=> b.ID }
        .flatMap { records ->
            check_unique_ids(records)
            records.collate(batch_size).withIndex()
        }
        .map { batch, index ->
            [ "${prefix}_batch_${String.format('%05d', index)}.json", groovy.json.JsonOutput.toJson(batch) ]
        }
}

def download_route(meta) {
    if (meta.Platform == "ILLUMINA") return "illumina_to_unpack"
    if (meta.Platform == "ONT" && meta.ont_format == "fastq") return "ont_format_fastq"
    if (meta.Platform == "ONT" && meta.ont_format in ["pod5", "fast5"]) return "ont_format_unbasecalled"
    return "other"
}

workflow IRODS_QUERY {
        main:
        BEEFEATER_SEARCH()
        | flatten() //beefeater emits one file per platform, so emit them one at a time
        | splitJson() //split that sample row into metadata
        | set { meta_with_paths_ch }

        if (params.save_metadata) {
            meta_with_paths_ch
            | collectFile() { map -> [ "lane_metadata.jsonl", groovy.json.JsonOutput.toJson(map) + '\n' ] }  // json, not map.toString(), so values holding ',' '=' or a newline survive
            | set{ metadata_only }

            METADATA_QUERIED(metadata_only, "irods_queried")
        }

        emit:
        meta_with_paths_ch
}

workflow IRODS_GET {
    /*
    Download search records in batches so unpacking can start as soon as the first batch lands,
    rather than once every file has downloaded
    */

    take:
    records_ch // map: one beefeater search record

    main:
    records_ch.branch { meta ->
            ont: meta.Platform == "ONT"
            other: true
        }
        | set { records_by_platform }

    batch_records(records_by_platform.other, params.irods_batch_size, "reads")
    | mix(batch_records(records_by_platform.ont, params.irods_ont_batch_size, "ont"))
    | collectFile { name_and_json -> name_and_json }
    | BEEFEATER_GET

    BEEFEATER_GET.out.downloaded
    | flatMap { downloaded_json, downloads ->
        def files = [downloads].flatten()
        def records = new groovy.json.JsonSlurper().parseText(downloaded_json.text)

        // downloads finish in any order, so pair each record with its file by name
        records.collect { meta ->
            def name = meta.local_path.tokenize('/').last()
            def downloaded_file = files.find { it.name == name }
            if (!downloaded_file) {
                log.warn "${meta.ID}: ${name} is listed in ${downloaded_json.name} but was not found among the downloads; skipping"
            }
            [ meta, downloaded_file ]
        }
        .findAll { meta, downloaded_file -> downloaded_file }
    }
    | set { downloaded_ch }

    emit:
    downloaded_ch // tuple val(meta), path(cram or collection directory)
}

workflow CRAM_EXTRACT {

    take:
    meta_with_cram_ch // tuple val(meta), path(cram)

    main:
    COLLATE_FASTQ(meta_with_cram_ch)

    if (params.cleanup_intermediate_files_irods_extractor) {
        COLLATE_FASTQ.out.files_to_remove.flatten()
            .filter(Path)
            .map { file -> NextflowTool.safeDelete(file, workflow.workDir, log) }
    }

    emit:
    reads_ch = COLLATE_FASTQ.out.fastq_channel // tuple val(meta), path(forward_fastq), path(reverse_fastq)
}

workflow IRODS_EXTRACTOR {

    take:
    input_irods_ch // map

    main:

    if (!params.search) {
        // Only records that are used below are downloaded
        input_irods_ch
        | filter { meta -> download_route(meta) != "other" }
        | IRODS_GET

        IRODS_GET.out.downloaded_ch.branch { meta, downloaded_file ->
                illumina_to_unpack: download_route(meta) == "illumina_to_unpack"

                ont_format_fastq: download_route(meta) == "ont_format_fastq"

                ont_format_unbasecalled: download_route(meta) == "ont_format_unbasecalled"
            }
        | set { downloaded_objects }

        // Extract the Illumina fastqs from the CRAMs and publish them to the output directory
        CRAM_EXTRACT(downloaded_objects.illumina_to_unpack)

        // Publish the reads and/or squiggles to the output directory
        PUBLISH_FASTQ(downloaded_objects.ont_format_fastq)
        PUBLISH_UNBASECALLED(downloaded_objects.ont_format_unbasecalled)

        illumina_reads_ch = CRAM_EXTRACT.out.reads_ch // tuple val(meta), path(forward_fastq), path(reverse_fastq)
        ont_reads_ch = PUBLISH_FASTQ.out.path_channel // tuple val(meta), path(fastq)
        ont_unbasecalled_ch = PUBLISH_UNBASECALLED.out.path_channel // tuple val(meta), path(fast5/pod5)
    } else {
        log.info "Search only mode enabled, no reads will be extracted or published."
        illumina_reads_ch = Channel.empty()
        ont_reads_ch = Channel.empty()
        ont_unbasecalled_ch = Channel.empty()
    }

    emit:
    illumina_reads_ch // tuple val(meta), path(forward_fastq), path(reverse_fastq)
    ont_reads_ch // tuple val(meta), path(fastq)
    ont_unbasecalled_ch // tuple val(meta), path(fast5/pod5)
}
