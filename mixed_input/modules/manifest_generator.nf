process MANIFEST_GENERATOR {
    label "cpu_1"
    label "mem_1"
    label "time_30m"

    publishDir "${params.outdir}", mode: 'copy'

    input:
    val input_dir

    output:
    path "manifest/manifest.csv", emit: 'ch_manifest_from_dir'

    script:
    manifest_script = "${projectDir}/assorted-sub-workflows/mixed_input/bin/generate_manifest.py"
    """
    mkdir -p manifest
    ${manifest_script} -i ${input_dir} -o manifest -v ${params.fastq_validation} -d ${params.max_depth}
    """
}
