# metawrap_qc

[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)


This workflow performs read quality control on raw metagenomic sequencing data.

## Pipeline Summary
Reads are trimmed of adapters and regions of low quality with TrimGalore. Host contamination is detected using BMTagger, with reads optionally retained separately for further analysis by adding the flag `--publish_host_data` .

The pipeline also generates stats about raw, trimmed and cleaned reads for each sample and collates these into one output CSV.

### Inputs
- Raw metagenomic paired-end reads (FASTQs) per sample

### Outputs
- Filtered, trimmed reads per sample (FASTQs)
- QC statistics CSV for all samples
- Host reads per sample (FASTQs) - OPTIONAL


### Dependencies
All dependencies including generate_stats.py and filter_reads.py scripts are available in the run time container. The dockerfile for this and the scripts themselves can be found [in the repo for metawrap_qc_nextflow](https://gitlab.internal.sanger.ac.uk/sanger-pathogens/pipelines/metawrap_qc_nextflow).


