# Preprocessing workflow

This workflow takes a pair of FASTQ files and performs a series of preprocessing steps to clean and prepare reads for any downstream analysis. Steps can include trimming, tandem repeat filtering (TRF), and host read removal, depending on the parameters yopu specify. By default trimming is turned on and TRF and host read removal is off. 

### Trimming

Trimming is turned on my default and conducted using Trimmomatic for adapter removal and quality trimming. Trimming ensures that only high-quality, informative sequences are used for downstream analysis. Note that adapter trimming is currently not carried out automatically in IROD by NPG. 

### Tandem repeat filtering (TRF)

If enabled by specifying `--run_trf true` Tandem Repeat Finder (TRF) is run to identify and remove reads with excessive tandem repeats (two or more adjacent, approximate copies of a pattern of nucleotides). This helps reduce artifacts caused by repetitive sequences that may interfere with mapping or assembly.

### Host read removal

When enabled by specifying `--run_bmtagger true`, BMTagger filters out reads that match the specified host genome (default is T2T-CHM13v2.0 human genome assembly). 

### Parameters
```
 General Preprocessing
      --publish_clean_reads
            default: true
            Save the pre-processed reads (gzip-compressed) in the preprocessing/ output folder
      --publish_trimmomatic_reads
            default: false
            Enable publishing of intermediate reads from trimmomatic process during pre-processing. Read sets will be uncompressed fastq files
      --publish_trf_reads
            default: false
            Enable publishing of intermediate reads from trf process during pre-processing. Read sets will be uncompressed fastq files
-----------------------------------------------------------------
 Trimmomatic
      --run_trimmomatic
            default: true
            Run Trimmomatic for adapter removal and trimming for quality
      --adapter_fasta
            default: /data/pam/software/trimmomatic/adapter_fastas/solexa-with-nextseqPR-adapters.fasta
            Path to fasta file containing adapter sequences
      --trim_window_size
            default: 4
            Sliding window size for read trimming
      --trim_baseq
            default: 20
            Average base quality cutoff for the sliding window for trimmomatic"
      --trim_min_length
            default: 70
            Minimum read length retained following trimming
      --trimmomatic_options
            default: ILLUMINACLIP:${params.adapter_fasta}:2:10:7:1 CROP:151 SLIDINGWINDOW:${params.trim_window_size}:${params.trim_baseq} MINLEN:${params.trim_min_length}
            Trimmomatic command line options
-----------------------------------------------------------------
 Tandem Repeat Finder (TRF)
      --run_trf
            default: false
            Run TRF for finding and removing tandem repeat from short reads
      --trf_cli_options
            default: 2 7 7 80 10 50 500 -h -ngs
            TRF command line options
-----------------------------------------------------------------
 BMTagger (Host Read Removal)
      --run_bmtagger
            default: false
            Run bmtagger for host read removal
      --publish_host_data
            default: false
            Publish the reads determined to originate from the host organism
      --bmtagger_db
            default: /data/pam/software/bmtagger
            Path to directory containing BMTagger database
      --bmtagger_host
            default: T2T-CHM13v2.0
            Reference genome version used for host read filtering
```