process BINETTE {
    tag "${meta.ID}"
    label 'cpu_8'
    label 'mem_64'
    label 'time_12'

    container 'quay.io/biocontainers/binette:1.2.1--pyh106432d_1'

    publishDir mode: 'copy', path: "${params.outdir}/${meta.ID}/binette/", pattern: "${report_txt}"

    input:
    // Four bin dirs in MAG_BINNING's fixed order. They are staged generically and
    // linked under their real binner's name below, so Binette's report says
    // "concoct_bins" when --no_gpu put CONCOCT in the first slot, never a
    // mislabelled "comebin_bins".
    tuple val(meta), path(bin_dirs, stageAs: 'bin_input_?'), path(assembly)

    output:
    // Emit the DIRECTORY, not a glob of its contents. `path("final_bins/*")`
    // yields a java.util.ArrayList of files, but the entire downstream metaWRAP
    // reassembly subworkflow assumes a single directory Path in THREE places:
    //   1. filesFromDir(dirPath)          -> .toFile().listFiles()
    //   2. COMBINE_BINS script            -> `for i in $(ls ${bins}); do cat ${bins}/$i ...`
    //   3. the bin-name parsers           -> basename.split('.')[1] as Integer
    // With a list, (1) throws, and (2) silently mis-expands to
    // `cat <all bins> <last bin>/$i`, which fails exit 1 AFTER writing a
    // hugely duplicated merged assembly - and because errorStrategy ignores it,
    // COLLECT_BINS/CHECKM then run on un-reassembled bins and the pipeline
    // reports success. Emitting the directory fixes all three at once.
    tuple val(meta), path("final_bins"), path(report_txt), emit: results

    script:
    report_txt = "${meta.ID}_binette_quality_report.tsv"
    names = [params.no_gpu ? 'concoct' : 'comebin', 'semibin2', 'metacat', 'maxbin2']
    assert bin_dirs.size() == names.size() : "BINETTE expected ${names.size()} bin dirs, got ${bin_dirs.size()}"
    links = [bin_dirs, names].transpose().collect { d, n -> "ln -s ${d} ${n}_bins" }.join('\n    ')
    bin_args = names.collect { "${it}_bins" }.join(' ')
    """
    ${links}

    binette \\
        --bin_dirs ${bin_args} \\
        --fasta_extensions .fasta .fa .fna \\
        --contigs ${assembly} \\
        --checkm2_db ${params.checkm2_db} \\
        -o binette_out -t ${task.cpus}

    mv binette_out/final_bins final_bins
    mv binette_out/final_bins_quality_reports.tsv ${report_txt}

    # Rename bins to the `bin.N.fa` convention that the downstream metaWRAP
    # reassembly subworkflow parses. reassembly.nf does
    #     basename.split('\\.')[1] as Integer      // "bin.9.fasta" -> 9
    # on both bin fastas and split-read filenames. Binette's native
    # `binette_binN.fa` makes parts[1] == "fa", so with reassembly enabled
    # (which is the DEFAULT - skip_reassembly = false) the run dies at DAG
    # construction with:
    #     ERROR ~ For input string: "fa"
    # Renaming here is the least invasive fix: it keeps Binette untouched and
    # satisfies both of reassembly.nf's filename parsers.
    ( cd final_bins && for f in binette_bin*.fa; do
        [ -e "\$f" ] || continue
        n=\${f#binette_bin}
        mv "\$f" "bin.\${n}"
      done )
    # keep the published quality report's bin names consistent with the files
    sed -i 's/^binette_bin/bin./' ${report_txt}
    """
}
