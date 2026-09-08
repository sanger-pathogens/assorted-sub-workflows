process COMEBIN {
    tag "${meta.ID}"
    label 'gpu'
    label 'cpu_5'
    label 'mem_8'
    label 'time_12'

    container 'quay.io/sangerpathogens/cuda_comebin:1.1.0'

    input:
    tuple val(meta), path(bam), path(bai), path(assembly)

    output:
    tuple val(meta), path(bins_out), emit: bins

    script:
    bins_out = "comebin/comebin_res/comebin_res_bins"
    """
    mkdir -p bamfiles
    ln -s "\$PWD/${bam}" "bamfiles/${bam}"
    ln -s "\$PWD/${bai}" "bamfiles/${bai}"

    run_comebin.sh -a ${assembly} -p bamfiles -o comebin -t ${task.cpus} -b 256
    """
}
