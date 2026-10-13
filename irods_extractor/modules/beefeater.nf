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


process BEEFEATER_SEARCH {
    label 'cpu_1'
    label 'mem_2'
    label 'time_12'

    container 'quay.io/sangerpathogens/beefeater:v1.3.0'

    output:
    path("*_output_*.json"), emit: json_ch // one file per platform, e.g. <date>_output_ILLUMINA.json

    script:
    // Get the absolute path to the output directory
    pwd = file(System.getenv("PWD"))
    outDirAbsPath = pwd.resolve(params.outdir).normalize().toString()

    // If the directory is nested feed the top level out dir to beefeater, if preprocessing and flat then point to preprocessing
    def dir_structure = ("${params.save_method}" == "nested") ? "${outDirAbsPath}" : "${outDirAbsPath}/${params.preexisting_output_dir}"

    // Options required for preventing the redownload of Illumina reads
    def illumina_read_output_directory = params.prevent_redownloads ? "--illumina_read_output_directory ${dir_structure}" : ""
    def illumina_publish_output = "--illumina_publish_output ${params.preexisting_output_dir}"
    def illumina_id_suffix = "--illumina_id_suffix ${params.existing_output_id_suffix}"
    def illumina_output_extension = "--illumina_output_extension ${params.existing_output_extension}"
    def illumina_publish_structure = "--illumina_publish_structure ${params.save_method}"
    prevent_redownload_opts = params.prevent_redownloads ? "${illumina_read_output_directory} ${illumina_publish_output} ${illumina_id_suffix} ${illumina_output_extension} ${illumina_publish_structure}" : ""

    // Options required for the download of ONT reads
    def ont_options = ""
    if (params.read_type.toLowerCase() == "ont") {
        if (params.ont_get_read_type == "reads") {
            ont_options = "--pick_best_ont"
        } else if (params.ont_get_read_type == "unbasecalled") {
            ont_options = "--unbasecalled"
        }
    }

    def query    = generate_query(params)

    def manifest = ""
    if (params.manifest_of_lanes) {
        def orig_manifest  = file(params.manifest_of_lanes)
        def lines          = orig_manifest.readLines()
        def new_header     = lines[0].split(',').collect { translateKey(it.trim()) }.join(',')
        def translated_csv = file("translated_manifest_of_lanes.csv")
        translated_csv.text = ([new_header] + lines.drop(1)).join('\n') + '\n'
        manifest = "--manifest ${translated_csv}"
    }
    def service_user = params.service_user ? "-su " : ""
    """
    beefeater search \\
        --file_format json \\
        --pool_size ${params.irods_pool_size} \\
        ${service_user} \\
        ${query} \\
        ${manifest} \\
        ${prevent_redownload_opts} \\
        ${ont_options}

    """
}


process BEEFEATER_GET {
    tag "${batch_json.baseName}"
    label 'cpu_1'
    label 'mem_2'
    label 'time_queue_from_normal'
    maxForks params.irods_max_concurrent_gets

    container 'quay.io/sangerpathogens/beefeater:v1.3.0'

    input:
    path(batch_json)

    output:
    tuple path("downloaded_*.json"), path("results/*"), emit: downloaded

    script:
    def service_user = params.service_user ? "-su " : ""
    """
    beefeater -v get \\
        --input ${batch_json} \\
        --file_format json \\
        --output_file downloaded.json \\
        --get_directory results \\
        --pool_size ${params.irods_pool_size} \\
        ${service_user}
    """
}
