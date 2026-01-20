process SPLIT_DEPTHS {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_100M'
    label 'time_30m'

    container 'quay.io/sangerpathogens/python-curl:3.11'

    input:
    tuple val(meta), path(depth_text)

    output:
    tuple val(meta), path(depth_out), emit: depths

    script:
    command = "${projectDir}/assorted-sub-workflows/mags_maker/binning/bin/split_depths.py"
    depth_out = "${meta.ID}_split_depths"
    """
    ${command} ${depth_text} ${depth_out}
    """
}

process MAXBIN2 {
    label 'cpu_2'
    label 'mem_1'
    label 'time_30m'

    container 'quay.io/biocontainers/maxbin2:2.2.6--h14c3975_0'

    input:
    tuple val(meta), path(depth_dir), path(assembly)

    output:
    tuple val(meta), path("maxbin2/"),  emit: bins

    script:
    """
    mkdir maxbin2

    run_MaxBin.pl -contig ${assembly} \\
        -markerset ${params.maxbin_markers} \\
        -thread ${task.cpus} \\
        -min_contig_length ${params.maxbin2_min_contig} \\
	    -out maxbin2/${meta.ID} \\
	    -abund_list ${depth_dir}/mb2_abund_list.txt

    #move stuff out of the bin that isn't to use
    #commented line is not added until later version

    mkdir maxbin_misc

    expected_files=(
        "${meta.ID}.marker"
        "${meta.ID}.noclass"
        "${meta.ID}.tooshort"
        "${meta.ID}.log"
        "${meta.ID}.summary"
        "${meta.ID}.seed"
        #"${meta.ID}.marker_of_each_bin.tar.gz"
    )

    for file in "\${expected_files[@]}"; do
        if [ -f "maxbin2/\${file}" ]; then
            mv "maxbin2/\${file}" maxbin_misc
        fi
    done
    """
}