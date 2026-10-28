# qc_mags

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.04.0-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

This workflow performs assembly quality control, with a focus on metagenome-assembled genomes (MAGs). Isolate genomes can be QC'd in the same run by declaring their type (see [Genome types: MAGs and isolates](#genome-types-mags-and-isolates)).

## Pipeline Summary

It runs the following steps on the assembly (`fasta` files):

1. **Input check**: every input FASTA is checked once before any QC. GenBank-style `//` lines and blank lines are removed, and sequences containing non-nucleotide characters are reported. See [Input checks](#input-checks).
2. **QUAST, CheckM2 and GUNC**: [QUAST](https://github.com/ablab/quast) evaluates genome/metagenome assemblies by computing various metrics. [CheckM2](https://github.com/chklovski/CheckM2) is used to predict the completeness and contamination of genomic bins (independent of their taxonomic classification). [GUNC](https://github.com/grp-bork/gunc) is also run to check for chimerism and contamination. This stage is the 'preqc'.
3. **MDMCleaner and seqkit**: [MDMCleaner](https://github.com/KIT-IBG-5/mdmcleaner) is used to remove contaminant contigs and provide reliable contig classification of metagenome assembled genomes (MAGs); isolates skip this step. Sequences that are less than a minimum length (see param `min_contig`) are then removed from every genome using [seqkit](https://bioinf.shenwei.me/seqkit/). See [Decontamination](#decontamination).
4. **QUAST, CheckM2 and GUNC**: The same commands as in stage 2 are run on the cleaned fasta outputs (postqc).
5. **gtdbtk classify_wf**: The cleaned assemblies are classified using steps detailed in the [gtdbtk docs](https://ecogenomics.github.io/GTDBTk/commands/classify_wf.html). See [Taxonomic classification](#taxonomic-classification).
6. **Reporting**: A summary TSV report per batch combines the QUAST, CheckM2 and GUNC results from before and after decontamination with the GTDB-Tk classification and a quality tier for each genome. All batch reports are also combined into one report for the run, with the manifest sample ID of each genome.

By default, each batch's intermediate work directories are deleted once the batch has finished. See [Work directory cleanup](#work-directory-cleanup).

### Parameters

Several parameters set paths to various databases used by pipeline tools:

- `--checkm2_db` (default: `/data/pam/software/checkm2_db/uniref100.KO.1.dmnd`)
- `--gunc_db` (default: `/data/pam/software/gunc/GTDB/gunc_db_gtdb95.dmnd`)
- `--mdmcleaner_db` (default: `/data/pam/software/mdmcleaner/v1`)
- `--gtdbtk_db` (default: `/data/pam/software/GTDBTk/release232_skani040`)

Other parameters:

- `--fasta_ext` (default: `fa`): Set the expected extension of fasta files in the input directories.
- `--genome_type` (default: `mag`): Genome type of every genome in the run, `mag` or `isolate`. See [Genome types: MAGs and isolates](#genome-types-mags-and-isolates).
- `--genome_types` (default: null): Optional TSV giving the type of individual genomes, for runs mixing MAGs and isolates. Genomes not listed use `--genome_type`. See [Genome types: MAGs and isolates](#genome-types-mags-and-isolates).
- `--fast_taxonomy` (default: false): Run GTDB-Tk without `--place_species`, using fewer CPUs and less memory. See [Taxonomic classification](#taxonomic-classification).
- `--report_config` (default: [`./assets/report_config.json`](./assets/report_config.json)). See [Report configuration](#report-configuration) for more details.
- `--autoqc_config` (default: null) Optional parameter to filter genomes that pass QC thresholds specified in the config. When applying this option you must either provide the path to a config or with the string "default" use the provided config at [`./assorted-sub-workflows/qc_mags/assets/autoqc_config.tsv`](./assorted-sub-workflows/qc_mags/assets/autoqc_config.tsv). See [AutoQC configuration](#autoqc-configuration) for more details.
- `--min_contig` (default: 1000) Contigs below this length will be filtered out by seqkit.
- `--skip_raw_reports` (default: false) Skip publishing raw reports from various QC tools, keeping only final summary reports.
- `--publish_mdmcleaner` (default: false) Also publish the MDMCleaner output FASTAs (before `min_contig` filtering) to `mdmcleaner/`.
- `--skip_cleanup` (default: false) Keep each batch's intermediate work directories. By default they are deleted once the batch has finished to save disk space, so a `-resume` reruns finished batches. See [Work directory cleanup](#work-directory-cleanup).
- `--batch_size` (default: 0) Number of genomes per batch. `0` processes each input element (sample) as its own batch. See [Batching large runs](#batching-large-runs).
- `--temp_file_storage` (default: `/tmp`): Specify a directory where GTDB-Tk can store temporary files during processing. Options are '/tmp', '/dev/shm' or 'null' (write to memory). See [GTDB-Tk runtime](#gtdb-tk-runtime) for more details.
- `--temp_space` (default: `30GB`): Request a specific amount of temporary working space to reserve for GTDB-Tk and MDMCleaner. See [GTDB-Tk runtime](#gtdb-tk-runtime) for more details.

## Configuration

### Report configuration

A configuration file can be supplied to the `--report_config` option to customise the final (per-batch) report generated by the pipeline. In particular, it allows (per-tool) customization of the identity column (to join the tool report tables together) and extract particular columns from these reports. The following tool names are currently mandatory: `QUAST`, `GTDBTK`, `GUNC`, `CHECKM2` (NOTE: tool name must be uppercase). In the following example config block, we select 2 columns (`classification` and `closest_genome_reference`) to extract from the GTDB-Tk report and to use `user_genome` as the identity column:

```
"GTDBTK": {
    "id_column": "user_genome",
    "keep_columns": [
        "classification",
        "closest_genome_reference"
    ]
}
```

Please refer to the default [`report_config.json`](./assets/report_config.json) for expected JSON structure.

NOTE: In addition to columns derived from the tool reports, every report includes the following columns:

| column                | description                                                                                                                                      |
| --------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------ |
| preqc_genome_name     | Name of the input fasta file minus extension                                                                                                     |
| postqc_genome_name    | Name of the fasta file output after processing with `seqkit` minus extension (`cleaned_<preqc_genome_name>_filtered_kept_contigs`)               |
| sample_or_strain_name | Name of the fasta file up to the last `_` character. Only reliable for `<sample>_<bin>` names; use `sample_id` in `all_final_report.tsv` instead |
| genome_status         | `mag` or `isolate`, from `--genome_type`/`--genome_types` (never guessed from the file name)                                                     |
| mimag_quality         | Quality tier `HQ`, `MQ` or `LQ` from post-QC CheckM2 and GUNC; empty if CheckM2 results are missing. See [Quality tiers](#quality-tiers)         |
| qc_review             | Set for isolates with post-QC completeness ≤ 90%, which may be complete genomes under-scored by CheckM2; empty otherwise                         |

`all_final_report.tsv` concatenates all batch reports and adds a `sample_id` column (second column), taken from the `ID` of the input element each genome came from. If the same genome file name appears under several IDs and the batch cannot tell them apart, `sample_id` lists all of them separated by `;`.

The default columns included in the final report are as follows:
_General_

- preqc_genome_name
- postqc_genome_name
- sample_or_strain_name
- genome_status
- mimag_quality
- qc_review

_CheckM2_

- checkm2_postqc_Completeness
- checkm2_postqc_Contamination
- checkm2_postqc_Genome_Size
- checkm2_preqc_Completeness
- checkm2_preqc_Contamination
- checkm2_preqc_Genome_Size

_GTDBTk_

- gtdbtk_classification
- gtdbtk_closest_genome_reference

_GUNC_

- gunc_postqc_n_contigs
- gunc_postqc_pass.GUNC
- gunc_preqc_n_contigs
- gunc_preqc_pass.GUNC

_QUAST_

- quast_preqc_\# contigs
- quast_preqc_Largest contig
- quast_preqc_Total length
- quast_preqc_GC (%)
- quast_preqc_N50
- quast_preqc_N90
- quast_preqc_auN
- quast_preqc_L50
- quast_preqc_L90
- quast*preqc*\# N's per 100 kbp
- quast*postqc*\# contigs
- quast_postqc_Largest contig
- quast_postqc_Total length
- quast_postqc_GC (%)
- quast_postqc_N50
- quast_postqc_N90
- quast_postqc_auN
- quast_postqc_L50
- quast_postqc_L90
- quast*postqc*\# N's per 100 kbp

A note on **number of contigs** metrics: QUAST "# contigs" reports all contigs (both in pre and post QC reports) whilst GUNC's "num_contigs" is defined as all contigs that contain a mapped gene.

Please see output reports for these individual tools to see which columns you can _add_ to your own configuration. Remove any by omitting them from your custom config. If a custom config leaves out `Completeness` or `Contamination` for CheckM2, `mimag_quality` is left empty; if it leaves out `pass.GUNC` for GUNC, no MAG can be `HQ` (a warning is logged in both cases).

If no genome in a batch has a GTDB-Tk classification, the report step fails rather than writing a report without taxonomy, since this means the GTDB-Tk summary was empty or did not match the genome names. Genomes missing a classification individually are listed as a warning in the task log.

If using `autoqc` ensure your report config includes the column you wish to filter on. See below for more information.

### Quality tiers

`mimag_quality` uses post-QC CheckM2 completeness and contamination with the same thresholds for every genome type:

| Tier | CheckM2 criteria                           |
| ---- | ------------------------------------------ |
| HQ   | completeness > 90% and contamination < 5%  |
| MQ   | completeness ≥ 50% and contamination < 10% |
| LQ   | anything else                              |

GUNC (post-QC) is then applied depending on the genome type:

- **MAGs:** a GUNC fail makes the genome `LQ`, and `HQ` additionally requires a GUNC pass, so a MAG without a GUNC result is `MQ` at most.
- **Isolates:** GUNC detects binning chimeras, which pure cultures cannot have, so it is reported but never changes the tier. Isolates at completeness ≤ 90% are flagged in `qc_review`, as CheckM2 is weakest on rare and novel species.

The rRNA/tRNA criteria of the full MIMAG high-quality standard are not assessed.

### AutoQC configuration

It is also possible to automatically filter genomes based on QC metrics produced in the summary report described in [Report configuration](#report-configuration).

To apply this function, you must enable AutoQC by supplying a configuration file path to `--autoqc_config` (config as described [here](../mixed_input/README.md#filter_metadatapy)) or by supplying the string "default" so that the default [autoqc_config.tsv](./assets/autoqc_config.tsv) is used.

The default AutoQC config's pass criteria are:

- Over 90% Completeness
- Under 5% Contamination
- Passes GUNC contamination checks

These criteria are applied to every genome, including isolates, so an isolate that fails GUNC also fails the default autoQC even though its `mimag_quality` can be `HQ`. To filter on the tier instead, use a custom config such as:

```tsv
column	filter	datatype
mimag_quality	"mimag_quality == 'HQ'"	str
```

It is possible to add filtering rules to the autoqc config for columns\* within the report config. If using a custom report config you can customise your AutoQC config: for the column using the tool name (lowercase) and `keep_column` from the report config seperated by an underscore. i.e. `<tool>_<keep_column>`. Note that for adding filters based on CheckM2 or GUNC you must specify whether to use preqc or postqc reports i.e. `<tool>_<qc_stage>_<keep_column>`, as shown in the filterable column list above.

\*NOTE: Unfortunately QUAST report columns are not currently columns you can AutoQC on. This is an area for future development and should be possible in future versions.

Please note - if you add the parameter `--autoqc_config` but do not supply either a valid config file path or "default" the pipeline will error on not being able to find a config file to use.

### Genome types: MAGs and isolates

The genome type decides how each genome is decontaminated and tiered:

|                       | MAG                                                  | Isolate                                                                                    |
| --------------------- | ---------------------------------------------------- | ------------------------------------------------------------------------------------------ |
| MDMCleaner            | Run: removes contaminant contigs (binning artefacts) | Skipped: cleaning a pure culture could only remove real accessory content such as plasmids |
| seqkit `--min_contig` | Run                                                  | Run                                                                                        |
| GUNC in the tier      | Gates the tier (see [Quality tiers](#quality-tiers)) | Reported only                                                                              |

Both types end up with the same output name, `cleaned_<genome>_filtered_kept_contigs.<ext>`.

The type must come from what you know about where each genome came from; it is never guessed from the file name. For a run of one type, set `--genome_type mag` (the default) or `--genome_type isolate`. For a run mixing both, also give `--genome_types`, a TSV with exactly this header and one row per genome:

```tsv
genome	genome_type
sample1_bin.3	mag
GCF_000005845.2	isolate
```

`genome` is the input file name without its extension (the report's `preqc_genome_name`), and `genome_type` is `mag` or `isolate` (case-insensitive). Genomes not listed use `--genome_type`. The table is checked when the workflow starts, which stops with an error naming the line for a wrong header, a row without exactly 2 tab-separated columns, an unknown type, or a genome listed twice.

### Inputs

- (Metagenomic) assemblies per sample. Channel emitting elements `[meta, [1.fa, 2.fa, ...]]`, where meta is a map containing metadata (including the mandatory `ID`, tracking a sample ID).

### Input checks

Before any QC, every input FASTA (plain or gzipped) is checked once per batch:

- GenBank-style `//` record terminator lines are removed
- blank (whitespace-only) lines are removed
- sequence lines containing characters that are not IUPAC nucleotide codes are reported, but not changed

Files that need lines removed are replaced by a cleaned copy with the same name; all other files are passed on unchanged. Every fix and warning is listed in `report/input_fasta_fixes.tsv` (columns: batch, file name, fix type, line count), and the number of fixes is logged as a warning at the end of the run.

### Outputs

Outputs are produced per batch: `<batch>` is the input element's `ID` when `--batch_size` is `0`, or `batch_1`, `batch_2`, ... when it is set.

- (Metagenomic) assemblies that have been decontaminated using MDMCleaner (MAGs only) and filtered for contig length according to `min_contig`, in `seqkit/` as `cleaned_<genome>_filtered_kept_contigs.<ext>`. If using `autoqc`, those passing are also published to `pass/fastas/`.
- MDMCleaner output FASTAs in `mdmcleaner/`, only with `--publish_mdmcleaner`.
- Reports generated by CheckM2, GUNC and QUAST in `pre_qc/` and `post_qc/`, and by GTDB-Tk (bacterial and archaeal summaries combined) in `post_qc/gtdbtk/` (unless using `--skip_raw_reports`).
- In `report/`:
  - `<batch>_final_report.tsv`: final report per batch (with customisable columns) combining individual reports from various pipeline tools.
  - `all_final_report.tsv`: all batch reports combined, with a `sample_id` column.
  - `<batch>_filter_summary.tsv`: summary of contigs filtered out per bin. Please note that the small contigs (below `min_contig`) are calculated first, the contigs MDMCleaner-failed contigs are only those removed despite their length being above `min_contig`.
  - `input_fasta_fixes.tsv`: fixes and warnings from the [input check](#input-checks).
  - `mdmcleaner_skipped.tsv`: MAGs that failed MDMCleaner and have no post-QC result (see [Decontamination](#decontamination)).

### Decontamination

MDMCleaner runs on each MAG separately, whatever the batch size. Contig names containing `:` (common in JGI/IMG genomes) are made safe for MDMCleaner and restored in its output.

If MDMCleaner fails on a MAG for a reason other than running out of memory or time, the task is ignored so the rest of the run continues, and that MAG has no post-QC results or row in the final report. These MAGs are listed in `report/mdmcleaner_skipped.tsv` (columns: batch, genome), which is always written (header only if none failed), and their number is logged as a warning. Genomes with no contigs left after `min_contig` filtering are likewise dropped from the post-QC steps.

### Taxonomic classification

GTDB-Tk `classify_wf` runs on the cleaned genomes, after decontamination and length filtering. Its summary combines the bacterial and archaeal results. GTDB-Tk first compares each genome with the GTDB species representatives by average nucleotide identity (skani); genomes without a close match are placed in the reference tree to classify them.

By default the workflow runs GTDB-Tk with `--place_species`, which also places the genomes already classified by skani into the tree. `--fast_taxonomy` drops this option and runs GTDB-Tk with fewer CPUs (2 instead of 8, or 4 instead of 64 for batches over 100 genomes) and 64 GB memory. Genomes that match a species representative are still classified to species by skani, but they are no longer placed in the tree, so the placement-related columns of the GTDB-Tk summary can differ.

### GTDB-Tk runtime

- Runtime and memory can be reduced for GTDB-Tk classficiation by specifying the `--temp_file_storage` option with either `/tmp` or `/dev/shm`. Writing these files to disk instead of keeping everything in memory reduces peak RAM usage by up to 89% and can improve runtime by up to 10%. This often allows the job to run with a smaller memory request, meaning it can start faster on cluster schedulers.
- If you know the size of your samples you can request a specific temp memory using `--temp_space <XX>GB`. This lets you reserve a specific amount of temporary memory/disk space. Typically you should request < 100 GB, and no more then 1000 GB, as larger requests may cause jobs to remain pending. Note that due to a know bug reported in the farm documentation, request half the memory you require as LSF double-accounts /tmp use see [here](https://ssg-confluence.internal.sanger.ac.uk/spaces/FARM/pages/101361225/Useful+LSF+resources#UsefulLSFresources-Resources).
- At most 10 GTDB-Tk and 10 MDMCleaner tasks run at once.

### Work directory cleanup

By default, once a batch is finished (its final report and contig filtering summary are written, plus its autoQC outputs if `--autoqc_config` is set), the work directories of that batch's tasks are deleted to save disk space. Only directories inside the run's `work/` directory are deleted; links to input files are removed without touching the files they point to.

Because the deleted tasks are no longer cached, `-resume` reruns any batch that had already finished. To keep everything, for example to debug a task or to resume a run cheaply, use `--skip_cleanup`.

### Dependencies

All software dependencies are containerised. Databases must be made available for CheckM2, GUNC, MDMCleaner and GTDBTk, paths supplied by params if not using the defaults.

### Troubleshooting

If a process in the pipeline fails with a non-zero exit code that is not recognised as a signal that memory or runtime limits were reached, the process will not be re-run and will be ignored. To understand why it has been ignored, please check the process-specific log file `.command.log` in the work directory of the affected process. With the default cleanup the work directories of finished batches are deleted, so rerun with `--skip_cleanup` to keep them.

Genomes missing from the final report are usually explained by `report/mdmcleaner_skipped.tsv` (MAGs that failed MDMCleaner) or `report/input_fasta_fixes.tsv` (problems found in the input FASTAs).

On the Sanger farm (HPC), a post-run script is available to automatically show the end of these log files (where informative messages can often be found) and return affected sample identifiers. See [here](https://ssg-confluence.internal.sanger.ac.uk/spaces/PaMI/pages/181078206/General+pipeline+info#Generalpipelineinfo-Usingthepipelinetracefile) for more details.

### Batching large runs

Many of the tools within the pipeline benefit from processing multiple MAGs simultaneously. To enable cross-sample batching, specify a `batch_size` using the flag:

```bash
--batch_size <INT>
```

This flag aggregates MAGs across samples into batches of N, with any remainder assigned to a final batch. For example, to define batches of 2000:

```bash
--batch_size 2000
```

Batch size involves a trade-off between throughput and resource usage. Larger batches improve computational efficiency but increase run-time per job, while smaller batches can reduce latency — particularly on an uncontested cluster — but at the cost of significantly higher memory overhead. Note that batches below ~1000 MAGs may incur disproportionate database loading overhead relative to actual QC compute time. A batch size of 1000–2000 MAGs is recommended for most use cases to strike a balance.

With batching, outputs are named `batch_1`, `batch_2`, ... instead of by sample. The input `ID` of each genome is recorded before batching, and is reported in the `sample_id` column of `report/all_final_report.tsv`.
