process MANIFEST_GENERATOR {
    label "cpu_1"
    label "mem_1"
    label "time_30m"

    container 'quay.io/sangerpathogens/manifest_generator:2.0.1-86d7316'

    publishDir "${params.outdir}/manifest", mode: 'copy'

    input:
    val(input_dir)

    output:
    path(output_manifest), emit: 'ch_manifest_from_dir'

    script:
    output_manifest = "manifest.csv"
    """
    manifest_generator.py --from_dir_recursive ${input_dir} --name ${output_manifest} --max_depth ${params.max_depth}
    """
}
