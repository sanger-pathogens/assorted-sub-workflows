process SOURMASH_SKETCH {
    tag "${meta.id}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_queue_from_normal'

    container 'quay.io/sangerpathogens/sourmash:4.5.0--hdfd78af_0'
    // more recent options are available:
    // container 'quay.io/sangerpathogens/sourmash:4.8.8--hdfd78af_0'

    input:
    tuple val(meta), path(merged_fastq)

    output:
    tuple val(meta), path(sourmash_sketch), emit: sketch

    script:
    sourmash_sketch="${meta.id}_sourmash_sketch"
    """
    sourmash sketch dna -p scaled=10000,k=31 ${merged_fastq} -o ${meta.id}_sourmash_sketch
    """
}

process SOURMASH_GATHER {
    tag "${meta.id}"
    label 'cpu_1'
    label 'mem_4'
    label 'time_queue_from_normal'

    container 'quay.io/sangerpathogens/sourmash:4.5.0--hdfd78af_0'
    // more recent options are available:
    // container 'quay.io/sangerpathogens/sourmash:4.8.8--hdfd78af_0'

    input:
    tuple val(meta), path(sourmash_sketch)

    output:
    tuple val(meta), path(sourmash_genomes), emit: sourmash_genomes

    script:
    sourmash_genomes="${meta.id}_sourmash_genomes.txt"
    """
    sourmash gather --dna ${sourmash_sketch} ${params.sourmash_db_abundance_estimation} -o sourmash.out
    # get species names out of sourmash output
    tail -n +2 sourmash.out | awk -F "," '{ print \$10 }' | sed 's|[][]||g' | sed 's|"||g' | awk '{ print \$1 }' > ${meta.id}_sourmash_genomes.txt
    """
}