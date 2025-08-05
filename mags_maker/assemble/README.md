# assemble

[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)


This script performs metagenomic binning of assembled contigs using MetaBAT2 (default), CONCOCT, and/or MaxBin2, leveraging multiple samples to improve bin quality.

## Workflow Summary

By default, this workflow assembles sequencing reads using metaSPAdes v3.14. Reads that do not map back to the resulting contigs are then re-assembled using MEGAHIT, which is better suited for recovering low-coverage regions. Users can choose to run only metaSPAdes or only MEGAHIT by specifying the --metaspades or --megahit parameters, respectively. The final output consists of the combined and sorted assemblies, with short contigs filtered out.

### Parameters
- `--assembly_stats` (Generate assembly statistics default: false)
- `--min_contig` (Minimum threshold for contig size (kb) default: 1500)
- `--fastspades` (Run SPAdes assembly module only default: null)
- `--lock_phred` (Set PHRED quality offset to 33 (SPAdes option) default: null)
- `--metaspades` (Run metaSPAdes assembly default: true)
- `--megahit` (Run MEGAHIT assembly default: true)
- `--output_transposed` (Output transposed quast reports default: false)

### Input
- Paired-end reads (FASTQs) 

### Outputs
- Sorted assemblies

### Dependencies
All dependencies are contained in runtime containers.

### Further information

If using a sra-lite format reads, or other formats not containing quality scores, the pipeline will attempt to automatically override the phred quality offset to 33.
If the pipeline is unable to automatically detect the lack of quality score before assembly steps it will fail the sample, set --lock_phred to provide a --phred_lock 33 to spades to address these failiures case. 

