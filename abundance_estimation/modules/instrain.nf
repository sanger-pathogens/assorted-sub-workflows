process INSTRAIN_PROFILE {
    tag "${meta.id}"
    label 'mem_4'
    label 'time_queue_from_small_slow2'

    container 'quay.io/sangerpathogens/instrain:1.9.0'
    
    if (params.instrain_full_output_abundance_estimation) { publishDir path: "${params.results_dir}/${meta.id}/instrain/", mode: 'copy', overwrite: true, pattern: "*_instrain_output" }
    publishDir "${params.results_dir}/${meta.id}/instrain/", mode: 'copy', overwrite: true, pattern: '*.tsv'
    
    input:
    tuple val(meta), path(sorted_bam), path(stb_file), path(genome_file)

    output:
    path("${meta.id}_instrain_output"), emit: full_output, optional: true
    tuple val(meta), path("${genome_info_file}"), emit: genome_info_file, optional: true
    path(sorted_bam), emit: sorted_bam
    tuple val(meta), path("${workdir}"), emit: meta_workdir

    script:
    genome_info_file="${meta.id}_genome_info.tsv"
    workdir="workdir.txt"
    """
    pwd > workdir.txt
    inStrain profile ${sorted_bam} ${genome_file} -o ${meta.id}_instrain_output -p ${task.cpus} -s ${stb_file} --skip_plot_generation ${params.instrain_profile_options} 2> instrain.err
    status=\${?}
    if [ \${status} -gt 0 ] ; then
        # try and catch known exceptions from the stored stderr stream
        grep "Exception: No paired reads detected" instrain.err && cat instrain.err 1>&2 && exit 7
        # if not caught known exception, process should not have exited yet - do it now spitting back the stored stderr and exit status
        cat instrain.err 1>&2 && exit \${status}
    else
        cat instrain.err 1>&2
    fi
    cp ${meta.id}_instrain_output/output/${meta.id}"_instrain_output_genome_info.tsv" ./${genome_info_file}
    """
}

process INSTRAIN_QUICKPROFILE {
    tag "${meta.id}"
    label 'mem_4'
    label 'time_queue_from_normal'

    container 'quay.io/sangerpathogens/instrain:1.9.0'
    
    if (params.instrain_quick_profile_abundance_estimation) { publishDir path: "${params.results_dir}/${meta.id}/instrain/", mode: 'copy', overwrite: true, pattern: "*_instrain_quick_profile_output" }
    publishDir "${params.results_dir}/${meta.id}/instrain/", mode: 'copy', overwrite: true, pattern: '*.tsv'
    
    input:
    tuple val(meta), path(sorted_bam), path(stb_file), path(genome_file)

    output:
    path("${meta.id}_instrain_quick_profile_output"), emit: quick_profile, optional: true
    path(sorted_bam), emit: sorted_bam
    tuple val(meta), path("${workdir}"), emit: meta_workdir

    script:
    workdir="workdir.txt"
    """
    pwd > "${workdir}"
    inStrain quick_profile ${sorted_bam} ${genome_file} -o ${meta.id}_instrain_quick_profile_output -p ${task.cpus} -s ${stb_file}
    """
}

process GENERATE_STB {
    label 'cpu_1'
    label 'mem_1'
    label 'time_queue_from_normal'

    input:
    tuple val(meta), path(sourmash_genomes)

    output:
    tuple val(meta), path("*.stb"), emit: stb_ch

    script:
    """
    sed 's|\$|_genomic.fna.gz|g' $sourmash_genomes > genomes.txt
    grep -w -f genomes.txt ${params.stb_file_abundance_estimation} > ${meta.id}_gtdb_subset.stb
    """
}

process GENERATE_INSTRAIN_SUMMARY {
    label 'cpu_1'
    label 'mem_1'
    label 'time_30m'

    publishDir "${params.results_dir}/abundance_summary", mode: 'copy', overwrite: true

    input:
    path(genome_info_files)

    output:
    path("instrain_summary*.tsv"),  emit: instrain_summary

    script:
    """
    ${projectDir}/bin/combine_instrain_output.sh "${params.custom_taxon_names_abundance_estimation}"
    """
}
