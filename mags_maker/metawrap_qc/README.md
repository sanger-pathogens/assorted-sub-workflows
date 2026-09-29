# metawrap_qc

[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)


This workflow performs read quality control on raw metagenomic sequencing data.

## Pipeline Summary
Reads are trimmed of adapters and regions of low quality with TrimGalore. Host contamination is detected using BMTagger, with reads optionally retained separately for further analysis by adding the flag `--publish_host_data` .

The pipeline also generates stats about raw, trimmed and cleaned reads for each sample and collates these into one output CSV.

### Parameters
- `--publish_host_data` (default: False)
- `--bmtagger_db` (default: "/data/pam/software/bmtagger")
- `--bmtagger_host` (default: "T2T-CHM13v2.0")

### Inputs
- Raw metagenomic paired-end reads (FASTQs) per sample

### Outputs
- Filtered, trimmed reads per sample (FASTQs)
- QC statistics CSV for all samples
- Host reads per sample (FASTQs) - OPTIONAL


### Dependencies

#### Scripts
All dependencies including generate_stats.py and filter_reads.py scripts are available in the run time container. The dockerfile for this and the scripts themselves can be found [in the repo for metawrap_qc_nextflow](https://gitlab.internal.sanger.ac.uk/sanger-pathogens/pipelines/metawrap_qc_nextflow).

#### Database
The `BMTAGGER` process requires access to the BMTagger database. If not using this pipeline on the Sanger HPC, you will need to biuld this database first.

Reference genome seqquences can be downloaded from the NCBI:
- for `hg38` please look for the latest annotated version here: https://hgdownload.soe.ucsc.edu/downloads.html
- for `T2T-CHM13v2.0`, use this NCBI RefSeq record: //ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz

To build the index, you will need the BMTagger executables; these are available from different sources:
- from source: [ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/bmtagger/]
- as a `conda` recipe from https://bioconda.github.io/recipes/bmtagger/README.html to make a `conda` environment
- using the same Docker image used in the pipeline: `quay.io/biocontainers/bmtagger:3.101--h470a237_4`; for this you need
  - Docker installed
  - to run the commands `bmtool` and `srprism` of the BMTagger tool suite having prepended this `docker run -v "$PWD:$PWD" -u "$(id -un):$(id -gn)" quay.io/biocontainers/bmtagger:3.101--h470a237_4 ...`

```sh
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/*fa.gz
# for T2T-CHM13v2.0, use `wget //ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz`
gunzip *fa.gz
cat *fa > hg38.fa
# for T2T-CHM13v2.0, use `cat *fa > T2T-CHM13v2.0.fa` and replace "hg38" with "T2T-CHM13v2.0" throughout
rm chr*.fa

bmtool -d hg38.fa -o hg38.bitmask
srprism mkindex -i hg38.fa -o hg38.srprism -M 100000 > index_hg38.o 2> index_hg38.e 
# this is a long process so can run this command in the background by appending a `&`
# the logs can then be watched for progress using `tail -f index_hg38.o`
```

