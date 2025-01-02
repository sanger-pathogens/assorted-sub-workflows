include { DOWNLOAD_METADATA; DOWNLOAD_FASTQS } from '../modules/ena_downloader'

workflow ENA_DOWNLOAD {
    take:
    accessions_file

    main:
    Channel.of([:])
    | combine(accessions_file)
    | set { download_metadata_input }

    DOWNLOAD_METADATA(download_metadata_input)
    | splitCsv(header: true, sep: "\t")
    | map { meta, full_metadata ->
        def sample_acc = full_metadata.sample_accession
        def cleaned_map = full_metadata.findAll { k, v -> v != '' }
        [ sample_acc, cleaned_map ] //staging sample_acc infront for groupTuple to output from ENADownloader
    }
    | set { sample_metadata }

    sample_metadata
    | map { id, metadata ->
        metadata.ID = id
        metadata
    }
    | filter { it.fastq_ftp.contains(';') } //if its paired its seperated by a semi-colon
    | map{ merged_meta ->
        def (read1_ftp, read2_ftp) = merged_meta.fastq_ftp.split(';')
        def read1_ftp_url = "ftp://${read1_ftp}"
        def read2_ftp_url = "ftp://${read2_ftp}"
        [ merged_meta, read1_ftp_url, read2_ftp_url ]
    }
    | DOWNLOAD_FASTQS
    | ifEmpty { error("Error: All Downloads failed") }
    | set { reads_from_ena_ch }

    emit:
    reads_from_ena_ch
}