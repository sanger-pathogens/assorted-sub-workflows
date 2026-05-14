process FASTQ2FASTA {
   label 'mem_1'
   label 'time_1'
   label 'cpu_1'
   tag {meta.ID}

   container 'quay.io/sangerpathogens/python-curl:3.11'
   /*
   *             Process: fastq_to_fasta

   This process converts FASTQ files into
   FASTA format using the `fastq2fasta.py` script.
   It takes as input FASTQ files along with metadata describing
   the sample. The output is the FASTQ in FASTA format.

   * --------------------------------------------------------------
   * Input:
      - tuple val(ID), val(read_type), val(meta), path(fastq)
      * ID: Unique identifier for the sample (e.g., sample ID).
      * read_type: Indicates whether the FASTQ file is for read 1, read 2, or unpaired reads (e.g., "1", "2", "unpaired")
      * meta: Metadata for the sample, which includes a
               unique identifier (e.g., sample ID) in `meta.id`.
      * fastq: Path to the FASTQ file.


   * Output:
      - tuple val(ID), val(read_type), val(meta), path(fastq), path(fasta)
      * ID: Unique identifier for the sample (e.g., sample ID).
      * read_type: Indicates whether the FASTQ file is for read 1, read 2, or unpaired reads (e.g., "1", "2", "unpaired")
      * meta: Same metadata as provided in the input.
      * fastq: Path to the original FASTQ file (same as input).
      * ${fasta}: FASTA file corresponding to the FASTQ file.

   * --------------------------------------------------------------
   * Dependencies:
      - Python3 required by `fastq2fasta.py` should be installed on
      the execution environment.

   * --------------------------------------------------------------
   */

   input:
   tuple val(ID), val(read_type), val(meta), path(fastq)

   output:
   tuple val(ID), val(read_type), val(meta), path(fastq), path(fasta), emit: converted

   script:
   fasta = "${meta.ID}_${read_type}.fasta"
   """
   ${params.script_src_path}fastq2fasta.py ${fastq} -o ${fasta}
    """
}