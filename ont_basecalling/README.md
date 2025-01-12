# Nanopore Basecalling

This subworkflow processes Nanopore sequencing data in POD5 or FAST5 formats, supporting file conversion, basecalling, demultiplexing, and quality control.

## Installation

1. [Install Nextflow](https://www.nextflow.io/docs/latest/install.html)

2. [Install Docker](https://docs.docker.com/engine/install/)

3. (Optional) Download the appropriate Dorado installer from the [repo](https://github.com/nanoporetech/dorado#installation). The path to the executable will be `<path to downloaded folder>/bin/dorado`

4. (Optional) Download the appropriate Dorado model from the [repo](https://github.com/nanoporetech/dorado/#available-basecalling-models)

   ```
   # Download all models
   dorado download --model all
   # Download particular model
   dorado download --model <model>
   ```

   If a pre-downloaded model path is not provided to the pipeline, the model specified by the `--basecall_model` parameter will be downloaded on the fly.

## Usage
Basic usage and mandatory parameters

```
nextflow run my-pipeline-importing-this-workflow/main.nf \
--raw_read_dir <directory containing FAST5/POD5 files> \
--additional_metadata <CSV mapping sample IDs to barcodes> \
--basecall_model <name of Dorado basecalling model> \
--barcode_kit_name <name of ONT barcoding kit used for multiplexing>
```

The CSV must be with the following stucture with each sample on a new line: 
ID,barcode,barcode_kit
21GUS-SR-008,01,SQK-NBD114-24

You can run the pipeline with `-profile laptop`. As well as enabling docker, the laptop profile allows the pipeline to be used offline by providing a local copy of a configuration file that is otherwise downloaded.

Should you need to run the pipeline offline, it is best to make use of pre-populated dependency caches. These can be created with any of the supported profiles (e.g. `-profile docker`) and involves running the pipeline once to completion. You will also need to provide a `--basecall_model_path` (see installation step 4).  
The config file can be edited to specify a local path for `--basecall_model_path` and `--dorado_local_path` in the laptop profile scope section, allowing quicker run configuration with `-profile laptop`.

You can override the default paths using the command line parameters directly when invoking nextflow or supplying an additional config file in which these parameters are set, using the `-c my_custom.config` nextflow option.

### Other parameters:

#### ONT Basecalling
Below are paramaters can be passed to customise the behaviour of the pipeline and their defaults
- --basecall = "true"
- --basecall_model = "dna_r10.4.1_e8.2_400bps_hac@v4.3.0"
- --basecall_model_path = ""
- --trim_adapters = "all"
- --barcode_kit_name = ["SQK-NBD114-24"] (currently this can only be edited via the config file)
- --read_format = "fastq"

#### Saving output files

- --save_fastqs = true