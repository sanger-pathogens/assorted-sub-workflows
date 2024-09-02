# kraken2bracken

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.04.0-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**kraken2bracken** is a pipeline for classifying bacterial reads using [Kraken2](https://github.com/DerrickWood/kraken2/) and obtaining relative abundance estimates with [Bracken](https://ccb.jhu.edu/software/bracken/index.shtml).

## Pipeline summary

**kraken2bracken** takes a Kraken2 database and sample manifest as input. It classifies reads in `.fastq` or `.fastq.gz` files with Kraken2. For each sample, it generates standard a Kraken2 report, sample level Kraken2 report and, optionally, files of classified and unclassified reads. As Kraken2 is run with `--report-minimizer-data` flag, the kraken reports will include additional columns (described [here](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#distinct-minimizer-count-information)). The pipeline then performs abundance estimation with Bracken, generating a standard Bracken report and a Kraken2 sample-format Bracken report for each sample.

If the given Kraken2 database does not contain a `database<read_len>mers.kmer_distrib` file and `--enable_building` is specified, it will generate one. If the file does not exist and `--enable_building` is not specified, the pipeline will fail with an error.

All relevant files are currently published in sample and process-specific directories within the supplied `--output` directory.

## Getting started

### Running on the farm (Sanger HPC clusters)

1. Load nextflow and singularity modules:
   ```bash
   module load nextflow ISG/singularity
   ```

2. Clone the repo:
   ```bash
   git clone --recurse-submodules git@gitlab.internal.sanger.ac.uk:sanger-pathogens/pipelines/kraken2bracken.git
   cd kraken2bracken
   ```

3. Start the pipeline:
   For example input, please see [Generating a manifest](#generating-a-manifest).

   Example:
   ```bash
   nextflow run . --input ./test_data/inputs/test_manifest.csv --read_len 150 --threshold 1 --classification_level 'G' --kraken2_threads 10 --outdir my_output
   ```

   It is good practice to submit a dedicated job for the nextflow master process (use the `oversubscribed` queue):
   ```bash
   bsub -o output.o -e error.e -q oversubscribed -R "select[mem>4000] rusage[mem=4000]" -M4000 nextflow run . --input ./test_data/inputs/test_manifest.csv --read_len 150 --threshold 1 --classification_level 'G' --kraken2_threads 10 --outdir my_output
   ```

   See [usage](#usage) for all available pipeline options.

4. Once your run has finished, check output in the `outdir` and clean up any intermediate files. To do this (assuming no other pipelines are running from the current working directory) run:

   ```bash
   rm -rf work .nextflow*
   ```

## Generating a manifest

Manifests supplied as an argument to `--input`, should be of of the following format:

```console
ID,R1,R2
test_id,./test_data/inputs/test_1.fastq.gz,./test_data/inputs/test_2.fastq.gz
```

Where column `ID` can be an arbitrary sample identifier, `R1` is a .fastq.gz file of forward reads, `R2` is the mate .fastq.gz file containing reverse reads.

Scripts have been developed to generate manifests appropriate for this pipeline:
- To generate a manifest from a file of lane identifiers visible to pf, use [this script](./scripts/generate_manifest_from_lanes.sh).
- To generate a manifest from a file of custom .fastq.gz paths, use [this script](./scripts/generate_manifest.sh).

Please run scripts with `-h` option for information on usage.

## Usage

```console
Usage:
    nextflow run main.nf
Options:
    --input                      Manifest containing per-sample paths to .fastq.gz files (mandatory)
    --kraken2_db                 Path to Kraken2 database (mandatory)
    --read_len                   Ideal length of reads in sample (mandatory)
    --kmer_len                   Length of kmers [Default: 35] (optional)
    --classification_level       Taxonomic rank to analyze for bracken2. Available options are 'D','P','C','O','F','G','S' [Default: 'S'] (optional)
    --threshold                  Minimum number of reads required for a classification at the specified classification_level [Default: 10] (optional)
    --get_classified_reads       Generate .fastq.gz files containing classified and unclassified reads for each sample [Default: False] (optional)
    --kraken2_threads            Threads to use for kraken2 [Default: 4] (optional)
    --bracken_threads            Threads to use for bracken [Default: 4] (optional)
    --outdir                     Specify output directory [Default: ./results] (optional)
    --help                       Print this help message (optional)
```

## Credits

kraken2bracken was produced by PAM informatics.

## Support

For further information or help, don't hesitate to get in touch via [path-help@sanger.ac.uk](mailto:path-help@sanger.ac.uk).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.
