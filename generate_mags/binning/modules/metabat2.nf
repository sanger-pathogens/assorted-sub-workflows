process CONTIG_DEPTHS {
    label 'cpu_2'
    label 'mem_100M'
    label 'time_30m'

    container 'quay.io/biocontainers/metabat2:2.12.1--1'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path(depth_text),  emit: depth

    script:
    depth_text = "${meta.ID}_depth.txt"
    """
    jgi_summarize_bam_contig_depths --outputDepth ${depth_text} ${bam}
    """
}

process CONTIG_DEPTHS_NO_INTRA {
    label 'cpu_1'
    label 'mem_100M'
    label 'time_30m'

    container 'quay.io/biocontainers/metabat2:2.12.1--1'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path(depth_text),  emit: depth

    script:
    depth_text = "${meta.ID}_depth.txt"
    """
    jgi_summarize_bam_contig_depths --outputDepth ${depth_text} --noIntraDepthVariance ${bam}
    """
}


process METABAT1 {
    label 'cpu_2'
    label 'mem_1'
    label 'time_30m'

    container 'quay.io/biocontainers/metabat2:2.12.1--1'

    input:
    tuple val(meta), path(depth_text), path(assembly)

    output:
    tuple val(meta), path("${meta.ID}_bin"),  emit: depth

    script:
    """
    metabat1 -i ${assembly} \\
        -a ${depth_text} \\
        -o ${meta.ID}_bin \\
        -m ${params.min_contig} \\
        -t ${task.cpus} \\
        --unbinned \\
        --seed ${params.bin_seed}
    """
}

process METABAT2 {
    label 'cpu_2'
    label 'mem_1'
    label 'time_30m'

    container 'quay.io/biocontainers/metabat2:2.12.1--1'

    input:
    tuple val(meta), path(depth_text), path(assembly)

    output:
    tuple val(meta), path("metabat/"),  emit: depth

    script:
    """
    metabat2 -i ${assembly} \\
        -a ${depth_text} \\
        -o metabat/${meta.ID}_bin \\
        -m ${params.min_contig} \\
        -t ${task.cpus} \\
        --unbinned \\
        --seed ${params.bin_seed}

    #move stuff out of the bin that isn't to use
    ## top 3 files are added in later version, we touch the files to avoid a crash 

    touch metabat/${meta.ID}_bin.BinInfo.txt . 
    mv metabat/${meta.ID}_bin.BinInfo.txt . 

    touch metabat/${meta.ID}_bin.lowDepth.fa .
    mv metabat/${meta.ID}_bin.lowDepth.fa .

    touch metabat/${meta.ID}_bin.tooShort.fa .
    mv metabat/${meta.ID}_bin.tooShort.fa .

    #Touch the files to avoid an error where no unbinned.fa file is made

    touch metabat/${meta.ID}_bin.unbinned.fa . 
    mv metabat/${meta.ID}_bin.unbinned.fa .

    # rename remaining fasta rather than fa
    for file in metabat/*.fa; do
        touch "\$file" "\${file%.fa}.fasta"
        mv "\$file" "\${file%.fa}.fasta"
    done
    """
}