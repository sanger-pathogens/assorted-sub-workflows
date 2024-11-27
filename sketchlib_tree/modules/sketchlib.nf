process SKETCH_ASSEMBLIES {
    label "cpu_2"
    label "mem_500M"
    label "time_30m"

    //will reach out to get a real fix rather than my attempt at fixing the index problem
    container 'quay.io/sangerpathogens/pp-sketchlib-rust:0.1.2_sd28_fix'

    input:
    tuple val(meta), path(assemblies)

    output:
    tuple val(meta), path("${sketch_db}.skm"), path("${sketch_db}.skd"), emit: assemblies_sketch

    script:
    sketch_db = "full_sketch"
    """
    sketchlib sketch -v -k 3,17,35 -o ${sketch_db} -s 1024 --seq-files ${assemblies}
    """
}

process SKETCH_CORE_ACC_DIST {
    label "cpu_2"
    label "mem_500M"
    label "time_30m"

    //will reach out to get a real fix rather than my attempt at fixing the index problem
    container 'quay.io/sangerpathogens/pp-sketchlib-rust:0.1.2_sd28_fix'

    input:
    tuple val(meta), path(skm), path(skd)

    output:
    tuple val(meta), path("total_core_data.tsv"), emit: distances

    script:
    """
    sketchlib dist -v ${skm} > total_core_data.tsv
    """
}
