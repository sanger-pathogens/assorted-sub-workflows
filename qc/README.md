# QC workflow

This workflow takes a pair of FASTQs and runs processes to determine the quality of a sample, outputting a `'pass'` or `'fail'` for each process.

### FastQC
If either of the FASTQs do not meet the required standard for the specified criteria in the FASTQC summary.txt file (based on the files provided to `--fastqc_pass_criteria` and `--fastqc_no_fail_criteria`), the sample is considered poor quality and a `'fail'` will be output. Otherwise, a `'pass'` will be output.

Example FastQC summary.txt:
```
PASS    Basic Statistics    SAMN11586388_1.fastq.gz
PASS    Per base sequence quality   SAMN11586388_1.fastq.gz
PASS    Per sequence quality scores SAMN11586388_1.fastq.gz
WARN    Per base sequence content   SAMN11586388_1.fastq.gz
WARN    Per sequence GC content SAMN11586388_1.fastq.gz
PASS    Per base N content  SAMN11586388_1.fastq.gz
PASS    Sequence Length Distribution    SAMN11586388_1.fastq.gz
PASS    Sequence Duplication Levels SAMN11586388_1.fastq.gz
PASS    Overrepresented sequences   SAMN11586388_1.fastq.gz
FAIL    Adapter Content SAMN11586388_1.fastq.gz
```

Example criteria file:
```
["Per base sequence quality", "Per sequence quality scores", "Per base N content"]
```

If provided to `--fastqc_pass_criteria`, these criteria would need a `'PASS'` for the sample to pass. If provided to `--fastqc_no_fail_criteria`, the criteria would need a `'PASS'` or a `'WARN'` for the sample to pass. Other criteria would be allowed to have a `'FAIL'` (or `'WARN'`).

### [kraken2bracken](../kraken2bracken/README.md)
If the abundance of the top (likely correct) genus or species in the Kraken style Bracken report is less than its specified threshold, the sample is considered contaminated/poor quality and a `'fail'` will be output. If the abundances meet the thresholds, a `'pass'` will be output.

### Parameters
```
 FastQC
      --save_fastqc
            default: false
            Flag to publish FastQC output
      --fastqc_pass_criteria
            default: path/to/repo/assorted-sub-workflows/qc/assets/fastqc_pass_criteria.json
            JSON file containing definition of an array specifying which items in the FastQC summary.txt are required to have the value PASS for the sample to be considered a pass
      --fastqc_no_fail_criteria
            default: path/to/repo/assorted-sub-workflows/qc/assets/fastqc_no_fail_criteria.json
            JSON file containing definition of an array specifying which items in the FastQC summary.txt are required to NOT have the value FAIL for the sample to be considered a pass (i.e. they could have WARN)

-----------------------------------------------------------------
 kraken2bracken
      --genus_abundance_threshold
            default: 90
            Fail the sample if the top genus abundance is lower than this

      --species_abundance_threshold
            default: 85
            Fail the sample if the top species abundance is lower than this
```