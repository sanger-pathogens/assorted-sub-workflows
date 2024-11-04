# QC workflow

This workflow takes a pair of FASTQs and runs processes to determine the quality of a sample, outputting a `'pass'` or `'fail'` for each process.

### FastQC
If there is a 'FAIL' in the FastQC summary.txt file for either of the FASTQs, the sample is considered poor quality and a `'fail'` will be output. If all stages have a 'PASS', a `'pass'` will be output.

### [kraken2bracken](../kraken2bracken/README.md)
If the abundance of the top (likely correct) genus or species in the Kraken style Bracken report is less than its specified threshold, the sample is considered contaminated/poor quality and a `'fail'` will be output. If the abundances meet the thresholds, a `'pass'` will be output.

### Parameters
```
 FastQC
      --save_fastqc
            default: false
            Flag to publish FastQC output

-----------------------------------------------------------------
 kraken2bracken
      --genus_abundance_threshold
            default: 90
            Fail the sample if the top genus abundance is lower than this

      --species_abundance_threshold
            default: 85
            Fail the sample if the top species abundance is lower than this
```