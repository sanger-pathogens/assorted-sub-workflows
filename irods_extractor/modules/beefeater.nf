def translateKey(in_key) {
    switch(in_key){ 
        case 'studyid':
        return 'study_id'

        case 'runid':
        return 'id_run'

        case 'laneid':
        return 'lane'

        case 'plexid':
        return 'tag_index'

        //if none of the cases return in_key
        default:
        return in_key
    }
}

def generate_query(params) {

    def keys = ['studyid', 'runid', 'laneid', 'plexid']

    def query_list = keys
        .findAll { params[it] != null }
        .collect { key ->
            def irods_key = translateKey(key)
            "-q ${irods_key}=${params[key]} "
        }

    return query_list.join('')
}


process BEEFEATER {
    label 'cpu_1'
    label 'mem_2'
    label 'time_1'

    container 'quay.io/sangerpathogens/beefeater:dev-v1.0.4'

    output:
    path("*output.csv"), emit: csv_ch //this is a json but the output file name is messed up to fix

    script:
    // Get the absolute path to the output directory
    pwd = file(System.getenv("PWD"))
    outDirAbsPath = pwd.resolve(params.outdir).normalize().toString()

    // If any preprocessing steps result in publishing _preprecessed.fastq.gz files then directory passed
    def preprocessing_params = [params.run_trimmomatic, params.run_trf, params.run_bmtagger]
    if (preprocessing_params.any { it == true } && params.prevent_redownloads_choice == "preprocessed") {
        final_directory = "preprocessing/"} else {
        final_directory = "${params.raw_reads_prefix}fastqs/"}
    // If the directory is nested feed the top level out dir to beefeater, if preprocessing and flat then point to preprocessing
    def dir_structure = ("${params.save_method}" == "nested") ? "${outDirAbsPath}" : "${outDirAbsPath}/${final_directory}"

    // Options required for preventing the redownload of Illumina reads
    def illumina_read_output_directory = params.prevent_redownloads ? "--illumina_read_output_directory ${dir_structure}" : ""
    def illumina_publish_output = "--illumina_publish_output ${final_directory}"
    def illumina_id_suffix = "--illumina_id_suffix ${params.existing_output_id_suffix}"
    def illumina_output_extension = "--illumina_output_extension ${params.existing_output_extension}"
    def illumina_publish_structure = ("${params.save_method}" == "nested") ? "--illumina_publish_structure nested" : "--illumina_publish_structure flat"
    prevent_redownload_opts = params.prevent_redownloads ? "${illumina_read_output_directory} ${illumina_publish_output} ${illumina_id_suffix} ${illumina_output_extension} ${illumina_publish_structure}" : ""


    def query    = generate_query(params)
    def manifest = params.manifest_of_lanes ? "--manifest ${params.manifest_of_lanes}" : ""
    def search   = params.search ? "" : "--get"

    """
    beefeater search \\
        --file_format json \\
        ${query} \\
        ${manifest} \\
        ${search} \\
        ${prevent_redownload_opts}
        
    """
}
