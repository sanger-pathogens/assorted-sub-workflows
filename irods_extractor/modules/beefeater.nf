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

    container 'quay.io/sangerpathogens/beefeater:v1.0.4'

    output:
    path("*output.csv"), emit: csv_ch //this is a json but the output file name is messed up to fix

    script:
    // Get the absolute path to the output directory
    pwd = file(System.getenv("PWD"))
    outDirAbsPath = pwd.resolve(params.outdir).normalize().toString()

    // If any preprocessing steps result in publishing _preprecessed.fastq.gz files then directory passed
    def preprocessing_params = [params.run_trimmomatic, params.run_trf, params.run_bmtagger]
    if (preprocessing_params.any { it == true }) {
        final_directory = "/preprocessing/"} else {
        final_directory = "/${params.raw_reads_prefix}fastqs/"}
    //If the directory is nested feed the top level out dir to beefeater, if preprocessing and flat then point to preprocessing
    def dir_structure = ("${params.save_method}" == "nested") ? "${outDirAbsPath}" : "${outDirAbsPath}/${final_directory}"
    def prevent_redownload = params.prevent_redownloads ? "--split_read_directory ${dir_structure}" : ""

    def query    = generate_query(params)
    def manifest = params.manifest_of_lanes ? "--manifest ${params.manifest_of_lanes}" : ""
    def search   = params.search ? "" : "--get"

    """
    beefeater search \\
        --file_format json \\
        ${query} \\
        ${manifest} \\
        ${search} \\
        ${prevent_redownload} 
        
    """
}
