process MANIFEST_GENERATOR {
    label "cpu_1"
    label "mem_1"
    label "time_10m"

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path input_dir

    output:
    path "manifest/manifest.csv", emit: 'ch_manifest_from_dir'

    script:
    """
    mkdir -p manifest
    ${projectDir}/scripts/generate_manifest.py \
        -i ${input_dir} \
        -o manifest \
        -v ${params.fastq_validation} \
        -d ${params.depth}
    """
}
