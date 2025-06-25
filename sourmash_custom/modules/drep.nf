process DREP_GENERATE_STB {
    label 'cpu_1'
    label 'mem_4'
    label 'time_queue_from_normal'
    publishDir "$params.outdir/Drep_stb", mode: "copy"

    container 'quay.io/biocontainers/drep:3.5.0--pyhdfd78af_0'

    input:
    path db_manifest

    output:
    path "${params.db_name}_drep.stb" , emit: stb_file

    script:
    """
    parse_stb.py --reverse -f ${db_manifest} -o ${params.db_name}_drep.stb
    """
}