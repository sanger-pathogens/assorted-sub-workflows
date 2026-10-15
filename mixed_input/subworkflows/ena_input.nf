include { DOWNLOAD_METADATA; DOWNLOAD_FASTQS } from '../modules/ena_downloader'
include { FILTER_EXISTING_OUTPUTS            } from '../../pipeline_chaining/subworkflows/skip_downloaded.nf'

workflow ENA_DOWNLOAD {
    take:
    accessions_file

    main:
    Channel.of([:])
    | combine(accessions_file)
    | set { download_metadata_input }

    DOWNLOAD_METADATA(download_metadata_input)
    | splitCsv(header: true, sep: "\t")
    | map { meta, full_metadata -> full_metadata.findAll { k, v -> v != '' } }
    | filter { row -> row.fastq_ftp?.contains(';') } //if its paired its seperated by a semi-colon
    | toList()
    | flatMap { rows ->
        // ENA gives one row per run. A sample with several runs would give them all the same ID and their
        // outputs would overwrite each other, so only then the run is added
        // (e.g. SAMD00162298_DRR171262, SAMD00162298_DRR171263); otherwise the ID is the sample accession
        def runs_per_sample = rows.countBy { row -> row.sample_accession }
        runs_per_sample.findAll { sample, n_runs -> n_runs > 1 }.each { sample, n_runs ->
            def run_ids = rows.findAll { row -> row.sample_accession == sample }.collect { row -> "${sample}_${row.run_accession}" }
            log.warn "ENA sample ${sample} has ${n_runs} paired runs; each is kept separately as ${run_ids.join(', ')}"
        }
        rows.collect { row ->
            def clash = runs_per_sample[row.sample_accession] > 1
            row + [ID: clash ? "${row.sample_accession}_${row.run_accession}".toString() : row.sample_accession]
        }
    }
    | map { merged_meta ->
        def (read1_ftp, read2_ftp) = merged_meta.fastq_ftp.split(';')
        def read1_ftp_url = "ftp://${read1_ftp}"
        def read2_ftp_url = "ftp://${read2_ftp}"
        [ merged_meta, read1_ftp_url, read2_ftp_url ]
    }
    | set { meta_readsftpurls_ch }
    
    if (params.only_new_input){
        FILTER_EXISTING_OUTPUTS(meta_readsftpurls_ch)
        FILTER_EXISTING_OUTPUTS.out.do_not_exist
        | set { meta_readsftpurls_to_process }
    } else {
        meta_readsftpurls_to_process = meta_readsftpurls_ch
    }

    DOWNLOAD_FASTQS(meta_readsftpurls_to_process)
    | set { reads_from_ena_ch }

    emit:
    reads_from_ena_ch
}