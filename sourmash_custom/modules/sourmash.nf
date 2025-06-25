process SOURMASH_SKETCH_DB { 
    label 'cpu_1'
    label 'mem_1' 
    label 'time_queue_from_normal'
    publishDir "$params.outdir/sourmash_db", mode: "copy"
    
    // container 'quay.io/biocontainers/sourmash:4.5.0--hdfd78af_0'
    container 'quay.io/sangerpathogens/sourmash:4.5.0--hdfd78af_0'

    input: 
    path(assembly_manifest)

    output:
    val(db_name), emit: db_name
    path(final_name)

    script:
    final_name = "${params.db_name}_s${params.sketch_size}k${params.klen}_sourmash.zip"
    """
    sourmash sketch dna -p scaled=${params.sketch_size},k=${params.klen} --from-file ${assembly_manifest} -o ${db_name}_s${params.sketch_size}k${params.klen}_sourmash.zip
    """

}