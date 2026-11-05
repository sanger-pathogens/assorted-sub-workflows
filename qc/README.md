# QC workflow

QC for paired short reads (`Platform` `ILLUMINA` or `PUBLIC_PAIRED`). Other platforms pass through unchanged.

1. **Human read removal** (bowtie2): matched to NPG's production method, using the NPG P4 container and NPG's
   T2T-CHM13v2.0 index. Each read is aligned single-end with `--local --very-sensitive-local`, and a pair is kept
   only if neither mate aligns. With `--host_removal auto` (default) it runs on every sample except reads from
   iRODS, which NPG has already depleted. This relies on `meta.input_source`, which `MIXED_INPUT` sets.
2. **fastp**: quality trimming matched to the previous Trimmomatic defaults (`SLIDINGWINDOW:4:20 CROP:151
   MINLEN:70`), plus forced polyG tail trimming.
3. **sylph**: each sample is sketched, then the sketches are profiled in batches (`--sylph_profile_batch_size`),
   so the database is loaded once per batch rather than once per sample. `sylph-tax` adds the taxonomy.
4. **MultiQC**: one report across all samples of the bowtie2 human alignment rates, fastp trimming and sylph-tax
   profiles (`--qc_multiqc_report`, default on). It is made once every sample has finished.

The fastp step replaces the [preprocessing](../preprocessing/README.md) subworkflow (Trimmomatic, BMTagger)
for pipelines using this QC. The previous QC workflow (FastQC, Kraken2/Bracken, pass/fail report) is kept in
[legacy_qc](../legacy_qc/README.md).

## Usage

```groovy
include { QC } from './assorted-sub-workflows/qc/qc.nf'

MIXED_INPUT()
| QC
```

and in `nextflow.config`:

```groovy
includeConfig "./assorted-sub-workflows/qc/qc.config"
```

## Outputs

```
<outdir>/
  preprocessed_fastqs/                  # QC'd reads (inside <ID>/ with --save_method nested)
    <ID>_preprocessed_{1,2}.fastq.gz
  <ID>/
    fastp/<ID>.fastp.{json,html}        # if --save_individual_fastp_reports
    sylph/<ID>_sylph_profile.tsv        # sylph profile (no file if sylph found nothing)
    sylph/<ID>.sylphmpa                 # sylph-tax MetaPhlAn-style profile
  preprocessing_summary_stats/
    fastp_statistics.tsv                # reads/bases before and after fastp, per sample
    host_removal_statistics.tsv         # reads in, aligned to human and kept, per sample (if any ran)
  multiqc/
    <date>-qc_report.html               # MultiQC report across all samples (if --qc_multiqc_report)
```

Emits: `preprocessed_reads_ch` (`[meta, read_1, read_2]`, plus any non-Illumina items unchanged), `fastp_stats`,
`host_removal_stats`, `sylph_report` and `sylphtax_mpa_report` (both `[meta, file]`), and `multiqc_report`.

Sylph profiling waits for every sample to be sketched before the batches start.

## Parameters

| Parameter                  | Default                              | Description                                                                                  |
| -------------------------- | ------------------------------------ | -------------------------------------------------------------------------------------------- |
| `--run_qc`                 | `true`                               | Run QC. When false, reads pass through unchanged and the QC outputs are empty.               |
| `--host_removal`           | `auto`                               | `auto` (all but iRODS reads), `all` or `none`.                                               |
| `--host_bowtie2_index`     | NPG T2T-CHM13v2.0                    | bowtie2 index of the human reference.                                                        |
| `--fastp_options`          | see `qc.config`                      | Options passed to fastp.                                                                     |
| `--fastp_trim_poly_g`      | `true`                               | Force polyG trimming (`-g`); false disables it (`-G`).                                       |
| `--save_individual_fastp_reports` | `true`                       | Publish fastp JSON/HTML reports.                                                             |
| `--publish_clean_reads`    | `true`                               | Publish QC'd reads.                                                                          |
| `--sylph_profile`          | `true`                               | Run sylph profiling.                                                                         |
| `--sylph_db`               | GTDB r232 c200 (`/data/pam/software`) | Sylph database.                                                                              |
| `--sylph_tax_metadata`     | GTDB r232                            | sylph-tax metadata matching `--sylph_db`.                                                    |
| `--sylph_profile_options`  | `-u`                                 | Extra `sylph profile` options; `-u` estimates the unknown fraction.                          |
| `--sylph_profile_batch_size` | `200`                              | Samples per `sylph profile` job.                                                             |
| `--save_sylph_sketches`    | `false`                              | Publish sylph sketches.                                                                      |
| `--qc_multiqc_report`      | `true`                               | Make one MultiQC report of bowtie2, fastp and sylph-tax results across all samples.          |
| `--qc_report_name`         | `qc_report`                          | MultiQC report name (`<date>-<qc_report_name>.html`).                                        |
| `--qc_multiqc_config`      | `qc/assets/multiqc_config.yaml`      | MultiQC config for the report (sample naming, section order).                                |

## Containers

| Software  | Version | Container                                          |
| --------- | ------- | -------------------------------------------------- |
| bowtie2 (NPG P4) | 2.5.4 | `ghcr.io/wtsi-npg/p4:0.48.5`                  |
| fastp     | 0.23.4  | `quay.io/biocontainers/fastp:0.23.4--hadf994f_2`    |
| sylph     | 1.0.0   | `quay.io/biocontainers/sylph:1.0.0--hb42e459_0`     |
| sylph-tax | 1.9.2   | `quay.io/biocontainers/sylph-tax:1.9.2--pyhdfd78af_0` |
