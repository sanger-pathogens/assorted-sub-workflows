# Assorted Sub-Workflows: a compendium of Nextflow modules for PaM Informatics bioinformatics pipelines

The following workflows are included, and can be integrated into any pipeline, in particular those maintained by the PaM Informatics team at the Wellcome Sanger Institute.

Usage and default parameters of each workflow are documented separately in their respective schema and config files.

## Combined short read input

The `COMBINED_INPUT` workflow enables input of short sequencing reads to your pipeline in multiple ways, which may be combined.


### Using read files already on disk

Input sequencing reads can be specified through a manifest file listing paths to pairs of fastq files. Manifests of read file paths can be supplied as an argument to `--manifest_of_reads` ; this manifest should be of the following format:
```
ID,R1,R2
genomeidA,/mydata/inputs/genomeidA_1.fastq.gz,/mydata/inputs/genomeidA_2.fastq.gz
genomeidB,/mydata/inputs/genomeidB_1.fastq.gz,/mydata/inputs/genomeidB_2.fastq.gz
```

Validation of the input manifest format occurs using the `INPUT_CHECK` subworkflow.

### Retrieving reads from iRODS

Sequencing data that was produced at the Wellcome Sanger Institute is natively stored on the iRODS platform, from where it can be readily obtained using the `IRODS_EXTRACTOR` subworkflow and its components `IRODS_QUERY` and `CRAM_EXTRACT`.  
Data can be queried from iRODS using all sorts of metadata fields that refer to specifics of the sequencing experiment and data to be retrieved.  
In practice, the user will usually rely on the core fields that define the "lane" id, or more correctly the lanelet id e.g. `48106_1#34`. These key fields are named (in the iRODS context): `study_id`, `id_run`, `lane` and `tag_index`; in the present context these properties were renamed for convenience as `studyid`, `runid`, `laneid` and `plexid`.  
Other search terms are supported however; see below.  

These search terms are used to build a query on iRODS that will select a dataset by restricting it to data objects (stored files) that match the parameter values. Each parameter restricts the set of data files that match and will be downloaded; when omitted, samples for all possible values of that parameter are retrieved.

#### Through the command-line arguments

The selected set of data files must be defined by a combination of the following CLI parameter options: `--studyid`, `--runid`, `--laneid`, `--plexid`, `--target` and `--type`, or through specification of the corresponding pipeline parameters through a Nextflow config file (passed with `-c` option).  

At least one of `studyid` or `runid` parameters must be specified. `laneid`, `plexid`, `target` and `type` are optional parameters that can be provided only in combination with `studyid` or `runid`; if these are specified without defining a `studyid` or `runid`, the request will be ignored (no iRODS data or metadata download) with a warning - this condition aims to avoid indiscriminate download of thousands of files across all possible runs.  

You can use a syntax like the below to specify a data set of reads present on iRODS (in this example, the workflow is being summoned by Nextflow pipeline called mypipeline):
```
mypipeline --studyid 7289 --runid 48106 --laneid 1 --plexid 34
```

#### Through a manifest of lane identifiers or other metadata terms
If you wish to retrieve multiple studies/runs at the same time, you can use the option `--manifest_of_lanes` to supply a manifest of lanes listing different combinations of the above iRODS metadata fields.

An example CSV with a few valid parameter combinations:
```
studyid,runid,laneid,plexid
7289,48106,1,36 
7289,48106,2,
,48106,,35
,48106,2,
```

Each row will result into a separate query to IRODS server.  
Note that some of these rows result in overlapping queries; while some files might be requested several times, the irods extractor workflow should only output one iteration of each file.  
Also, it is worth noting that iRODS querying can be slow, and that it is more efficient to do fewer queries that match several files instead of many queries, one for each read file set.

The order of parameters fields in the CSV is not relevant, but the user needs to specify at least one of the `studyid` or `runid` fields; leaving a parameter field blank means data matching all possible values will be selected. Again, `laneid` and `plexid` can only be used in combination with `studyid` or `runid`, or other fields - this is a rule we introduce to prevent accidental mass download that would ensue from not specifying the run.   
Similarly, one cannot submit an iRODS query solely based on `target` or `type` metadata tags, as this query would catch too many file objects.
For instance, these parameter combinations are invalid (run id not specified where lane or plex ids were) and will be ignored (no data or metadata download) with a warning:
```
studyid,runid,laneid,plexid,type,target
,,1,36,,
,,1,36,,0
,,,,bam,
```

Finally, using a manifest enables the saerch of a wider set of metadata fields, namely all metadata fields featured in the iRODS database system, as defined under the Dublin Core standard; see [here](https://github.com/wtsi-npg/irods-metadata/blob/master/irods_sample_metadata.md).

In practice, you can thus query samples based on attributes such as the organism they were collected from; this would be using the `sample_common_name` field, e.g.:  
```
sample_common_name,type
Romboutsia lituseburensis,cram
```
Note that, due to the large number of possible search terms used, no restriction rule is applied when using native iRODS search terms, so please be cautious when using this as it could result in very many file output (which might be what you desire!).

The workflow should issue a notice of how many files are considered for download, which may allow you to stop it in time if you get scared of the numbers.

### Saving (or not saving) Fastq files

When retrieving data from iRODS, one of the steps - enabled by the `IRODS_EXTRACTOR` subworkflow - is to download sequencing read files from iRODS (most often in CRAM format) and to convert them to the Fastq format. In most pipeline applications, the only purpose of these fastq files are to be passed on to further processing steps, and permanently saving fastq files on disk is not required. In fact, in most cases, it is not advisable as these files can be very large and will soon eat up all your disk space; this is why the workflow parameter `save_fastqs` defaults to `false`.  
It is the `irods_extractor` standalone pipeline use case, however, downlaoding fastqs is the purpose of the pipeline; see [the irods_extractor pipeline repo](https://gitlab.internal.sanger.ac.uk/sanger-pathogens/pipelines/irods_extractor).  
In that use case, or when using `IRODS_EXTRACTOR` subworkflow as part of another pipeline where you may want to retain the reads, you may consider the behaviour described below.

#### Updating the result folder for Fastq files

When attemting to download files from iRODS, the `IRODS_EXTRACTOR` subworkflow checks if the expected ouput fastq files exist in the `raw_fastq/` folder (default parameter value, folder name can be changed) under the main `results/` folder for each sample e.g. in `results/12345_1#67/raw_fastq/12345_1#67_1.fastq.gz`. If matching files are present, it then skips download and any further processing for this sample.  
The idea behind this mechanism is that when pulling reads from iRODS using a given set of query terms, the contents of the platform may evolve i.e. every so often new data files may become available that match the query, typically as new sequencing runs are completed as part of an ongoing study. To access those reads and generate any downstream analysis, there is no need to know the specifics of the new data to selectively download them. Instead, the whole study/dataset can be requested again and `IRODS_EXTRACTOR` will know to only update the results folder rather than re-downloading and re-processing everything - that is if you keep the same results folder and have the right values set for the relevant parameters: `preexisting_fastq_tag` and `split_sep_for_ID_from_fastq`, see help message for details; default behaviour is to update.  
Once it has processed the iRODS query and identified the set of files to be downloaded, the workflow will print the expected fastq file pair number to the screen (or standard output stream). It will then compare it the set of files already occurring under the `results/` folder and compute the difference, and will then print again how many files will actually be downloaded, which it will then proceed to do. You may have a look at these messages to make sure expected contents are downloaded. 
Note that in the case where you indeed desire to re-download the whole dataset and re-run the downstream analyses - for instance because a fix or new features were introduced in a new version of the pipeline - then you will want to turn this behaviour off, or to choose another destination for the results.

### Metadata search only


## Strain Mapper: short read mapping to a reference genome and variant genotype calling

The `STRAIN_MAPPER` workflow maps short read sequences to a given reference genome using a choice of mapping tools including **Bowtie2** and **BWA-MEM**. It also calls genotypes with **bcftools** (**GATK** will soon be supported as an alternative calling method) and thus generates a VCF file containing genotype  information, which is then used to create a consensus sequence (or "pseudo-genome") in Fasta format amalgamating variants to the reference sequence.

### Workflow details

#### Read mapping

Mapping can be done using Bowtie2 or BWA-MEM

#### Genotype calling

The final VCF file will include both REF and ALT alleles, or only variants if setting `only_report_alts = true`. By default, genotype calling is done using `bcftool mpileup` followed by `bcftools call`. A critical option is `minimum_base_quality` (default: 20), which will determine what bases are included at the pileup stage based on the basecalling quality (BQ); this has a significant impact in downstream genotype mapping quality and as a result, of how many sites are reported with a supported genotype, or represented with an `N` in the consensus sequence (see below).

Work in progress: inclusion of GATK as an alternative genotype caller is, enabling precise resolution of indels (see dev branch `PAT-1857_update_to_include_gatk_haplotypecaller`; MR [here](https://gitlab.internal.sanger.ac.uk/sanger-pathogens/pipelines/assorted-sub-workflows/-/merge_requests/14)).

#### Consensus sequence generation

a Python script is used to parse the VCF file and to report well-supported gentype calls into the context of the reference sequence, using `N` characters to represent sites where the genotype is uncertain due to too low mapping quality of read coverage.

#### Genotype/Variant filtering

By default, filtering using `bcftools view` removes genotype records (ALT and REF) with a mapping quality (MQ) score below 50 and ensures that they are represented on at least 3 forward and reverse reads. In addition there must be more than 8 reads in total supporting this variant call.

To edit these please follow the BCFTOOLS documentation as seen on the`bcftools view` page under "expressions": https://samtools.github.io/bcftools/bcftools.html#expressions


### Output 

The main output files are the following:

- the final VCF file (with filtered high-quality variants, or all variants if filtering as disabled), written into the folder `final_vcf/`
- the consesnsus sequence in Fasta format, written into the folder `curated_consensus/`

By default, most intermediate files are not published; when chosen to be published, they will be written in process-specific sub-folders under the sample-specific folder under the `results/` folder (or what was the supplied to the `--output` option). This includes:

- the raw VCF file from `bcftools call` before filtering low-quality variants. Use `keep_raw_vcf = true` to keep write it into `raw_vcf/` folder.
- the BAM alignment files produced by sorting with `samtools sort` and after deduplication with `picard MarkDuplicates` are discarded, unless setting `keep_sorted_bam = true` and `keep_dedup_bam = true` so that they are written to `samtools_sort/` and `picard/`, respectively.

In addition, optional analyses can be run and written alongside the above results file:
- `samtools stats` and `samtools flagstats` are run if `samtools_stats = true` (default) and results are saved in the `samtools_stats/` folder
- DeepTools `bamCoverage` is run to generate Bigwig coverage files, saved in the `bigwig/` folder
- The pipeline will automatically build reference and alignment indexes if it doesn't find them in the same directory as the supplied `--reference`, written into the folders `bowtie2/` or `bwa/` depending if the alignment of reads is made using Bowtie2 or BWA-MEM.

## Phylogenetic tree building