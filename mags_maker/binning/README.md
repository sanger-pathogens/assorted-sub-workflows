# assemble

[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)


This script bins contigs from an assembly using metaBAT2, CONCOCT and MaxBin2 into metagenomic bins.

## Workflow Summary

The script uses metaBAT2 (optional parameter --matabat1) CONCOCT and MaxBin2 to bin the contigs. The final output is sets of bins for each tool used.

### Parameters
- `--bin_seed` (Set a fixed seed for metaBAT1/2 and CONCOCT default: 1234)
- `--min_contig` (Minimum threshold for contig size (kb) default: 1500)
- `--maxbin_markers` (Marker set to choose from 107 marker genes present in >95% of bacteria, or 40 marker gene sets that are universal among bacteria and archaea default: 140??)
- `--metabat1` (Run metaBAT1 default: null, uses metaBAT2)
- `--keep_unbinned` (Keep unbinned contigs default: false)


### Input
- Sorted assemblies
- Paired-end reads (FASTQs) 

### Outputs
- Bin sets per sample per binning tool

### Dependencies
All dependencies are contained in runtime containers.

