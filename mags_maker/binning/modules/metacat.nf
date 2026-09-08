process METACAT {
    tag "${meta.ID}"
    label 'cpu_8'
    label 'mem_4'
    label 'time_12'

    container 'quay.io/sangerpathogens/metacat:1.0.6'

    input:
    tuple val(meta), path(bam), path(bai), path(assembly)

    output:
    tuple val(meta), path(bins_out), emit: bins

    script:
    bins_out = "metacat_out"
    """
    mkdir -p ${bins_out}
    MetaCAT coverage -b ${bam} -o ${bins_out}/coverage -tc ${task.cpus}
    MetaCAT seed -f ${assembly} -o ${bins_out}/seed -t ${task.cpus}
    MetaCAT cluster -f ${assembly} -c ${bins_out}/coverage -s ${bins_out}/seed -o ${bins_out}/metacat -t ${task.cpus}
    """
}
