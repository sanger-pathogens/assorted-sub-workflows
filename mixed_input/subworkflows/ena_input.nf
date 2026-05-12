include { DOWNLOAD_METADATA; DOWNLOAD_FASTQS } from '../modules/ena_downloader'
include { FILTER_EXISTING_OUTPUTS       } from './../../../pipeline_chaining/subworkflows/skip_downloaded.nf'

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
        full_metadata.ID = full_metadata.sample_accession
        def cleaned_metadata = full_metadata.findAll { k, v -> v != '' }
        cleaned_metadata
    }
    | filter { it.fastq_ftp.contains(';') } //if its paired its seperated by a semi-colon
    | map{ merged_meta ->
        def (read1_ftp, read2_ftp) = merged_meta.fastq_ftp.split(';')
        def read1_ftp_url = "ftp://${read1_ftp}"
        def read2_ftp_url = "ftp://${read2_ftp}"
        [ merged_meta, read1_ftp_url, read2_ftp_url ]
    }
    .set{ meta_readsftpurls_ch }
    
    if (params.only_new_input){
        FILTER_EXISTING_OUTPUTS(meta_readsftpurls_ch)
        FILTER_EXISTING_OUTPUTS.out.do_not_exist
        .set { meta_readsftpurls_to_process }
    } else{
        meta_readsftpurls_to_process = meta_readsftpurls_ch
    }

    DOWNLOAD_FASTQS(meta_readsftpurls_to_process)
    | ifEmpty { error("Error: All Downloads failed") }
    | set { reads_from_ena_ch }

    emit:
    reads_from_ena_ch
}