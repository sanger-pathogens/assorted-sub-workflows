params.script_src_path="${projectDir}/rvi_toolbox/bin/" 


process FASTQ2FASTA {
    /*
    *             Process: fastq_to_fasta

     This process converts paired-end FASTQ files into
     FASTA format using the `fastq2fasta.py` script.
     It takes as input two FASTQ files (e.g., forward
     and reverse reads) along with metadata describing
     the sample. The output consists of two corresponding
     FASTA files, one for each input FASTQ file.

    * --------------------------------------------------------------
    * Input:
       - tuple val(meta), path(fastq_1), path(fastq_2)
         * meta: Metadata for the sample, which includes a
                unique identifier (e.g., sample ID) in `meta.id`.
         * fastq_1: Path to the first FASTQ file (e.g., forward reads).
         * fastq_2: Path to the second FASTQ file (e.g., reverse reads).

    * Output:
       - tuple val(meta), path("${meta.id}_1.fasta"), path("${meta.id}_2.fasta")
         * meta: Same metadata as provided in the input.
         * ${meta.id}_1.fasta: FASTA file corresponding to
                              the first FASTQ file.
         * ${meta.id}_2.fasta: FASTA file corresponding to
                              the second FASTQ file.

    * --------------------------------------------------------------
    * Dependencies:
       - Python3 required by `fastq2fasta.py` should be installed on
       the execution environment.

    * --------------------------------------------------------------
    */
    tag {meta.id}
    container "quay.io/gsu-pipelines/rvi-vp-basecontainer"

    input:
        tuple val(meta), path(fastq_1), path(fastq_2) 

    output:
        tuple val(meta), path("${meta.id}_1.fasta"), path("${meta.id}_2.fasta")

    script:
    """
    ${params.script_src_path}fastq2fasta.py ${fastq_1} -o ${meta.id}_1.fasta
    ${params.script_src_path}fastq2fasta.py ${fastq_2} -o ${meta.id}_2.fasta
    """
}