# bin_refinement

[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)


This workflow refines metagenomic bins output by multiple binning tools (e.g. MetaBAT2, CONCOCT, MaxBin2) based on completeness and contamination metrics for each bin and all combinations.

## Pipeline Summary
All combinations of input bins (e.g. 2-wise, 3-wise...n-wise) are explored for potential merges. These combinations are refined by the tool Bin_refiner and quality of all bins - both input and refined - are assessed according to the minimum completeness and maximum contamination, by default with CheckM2 but optionally CheckM by using `--checkm1`. High quality bins are then merged and contigs are dereplicated, with ambiguous contigs assigned to the highest quality bin by default. Alternatively, ambiguous contigs can be removed entirely by setting (`--strict_bins`).


### Parameters
- `--checkm1` (default: False, uses CheckM2)
- `--checkm2_db` (default: "/data/pam/software/checkm2_db/uniref100.KO.1.dmnd")
- `--min_completeness` (default: 50)
- `--max_contamination` (default: 5)
- `--strict_bins` (default: False)

### Inputs
- Sets of bins per sample (output of binning tools)
- Paired-end reads (FASTQs) per sample??

### Outputs
- High-quality, dereplicated, refined and non-empty bins per sample
- CheckM/CheckM2 quality reports


### Dependencies
All dependencies are contained in runtime containers.


