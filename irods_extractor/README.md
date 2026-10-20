# iRODS extractor subworkflow

Retrieves sequencing data from IRODS system with
[beefeater](https://gitlab.internal.sanger.ac.uk/sanger-pathogens/beefeater), and turns it into reads:
Illumina CRAMs/BAMs are unpacked into paired `.fastq.gz` files, and ONT reads or squiggles are published as they are.

Pipelines normally reach this subworkflow through [`MIXED_INPUT`](../mixed_input/README.md), which calls
`IRODS_QUERY | IRODS_EXTRACTOR` whenever iRODS input is given. How to choose data can be done via many methods such as the `--studyid`, `--runid`,
`--laneid`, `--plexid` options, `--manifest_of_lanes`, and searching on other metadata fields

:warning: **iRODS access required.** You need an iRODS account and an active session (`iinit`) before running.
If you get a server error or a missing environment file message, request an account via
[https://sanger.freshservice.com/support/catalog/items/56](https://sanger.freshservice.com/support/catalog/items/56).

[[_TOC_]]

## How it works

```
IRODS_QUERY       BEEFEATER_SEARCH (1 job)  ->  one record per data object, with all its metadata
                  METADATA_QUERIED          ->  metadata_irods_queried_<datetime>.csv
IRODS_EXTRACTOR   route each record by data type; warn about any that cannot be handled
  IRODS_GET       split records into batches -> BEEFEATER_GET (one job per batch, a few at a time)
                  pair each downloaded file with its record
  CRAM_EXTRACT    COLLATE_FASTQ: Illumina CRAM/BAM -> <ID>_raw_{1,2}.fastq.gz
                  PUBLISH_FASTQ / PUBLISH_UNBASECALLED: ONT reads / squiggles
```

### 1. Search (`IRODS_QUERY`)

`BEEFEATER_SEARCH` runs `beefeater search` once, without downloading anything. It writes one JSON file per
platform, which is split into one record (a map of iRODS metadata plus `ID`, `Platform` and `irods_path`) per
data object. With `--save_metadata true` (default), all records are also saved to
`<outdir>/metadata_irods_queried_<datetime>.csv`.

With `--search true` the subworkflow stops here: nothing is downloaded and the emitted (output) read channels are empty.

> Narrow queries are slow on the iRODS side. Adding `lane`/`tag_index` to a query (CLI or manifest) took
> ~3 min for 3 plexes, against ~2 s for a whole study (study 6495, October 2026). Prefer `studyid`/`runid` where
> you can.

### 2. Routing by data type (`IRODS_EXTRACTOR`)

Each record is routed on its `Platform` (and in addition, for ONT, its `ont_format`):

| Record                                        | Route                                                                        |
| --------------------------------------------- | ---------------------------------------------------------------------------- |
| `Platform` `ILLUMINA`                         | downloaded, then unpacked to paired fastq (`CRAM_EXTRACT`)                   |
| `Platform` `ONT`, `ont_format` `fastq`        | downloaded, fastq files concatenated and published (`PUBLISH_FASTQ`)         |
| `Platform` `ONT`, `ont_format` `pod5`/`fast5` | downloaded and published as squiggles (`PUBLISH_UNBASECALLED`)               |
| anything else (e.g. PacBio, ONT `bam`)        | **not downloaded**; one warning per run lists the data types and example IDs |

Records that cannot be handled do not fail the run. To retrieve them, use the latest beefeater directly
(`beefeater search ... --get`).

### 3. Batched download (`IRODS_GET`)

Records are downloaded in batches, so unpacking starts as soon as the first batch has landed rather than after
every file has downloaded.

- Records are sorted by `ID` before batching, so the same records always give the same batches and `-resume`
  reuses completed downloads.
- Records that share an `ID` stop the run with an error listing them, as their downloads would overwrite each
  other.
- **Illumina (and other non-ONT) records** are split evenly into
  `max(1, min(irods_batch_count, n_records / irods_min_batch_size))` batches. With the defaults (20, 20):
  fewer than 40 records is one batch, 400 records or more is 20 batches, and beyond that batches grow rather than
  adding jobs.
- **ONT records** are collections that can be very large, so they are batched by a fixed number of collections
  per job (`irods_ont_batch_size`, default 1).
- Each batch is written to a JSON file and downloaded by one `BEEFEATER_GET` job running `beefeater get`. At most
  `irods_max_concurrent_gets` jobs run at once.
- Downloads finish in any order, so each downloaded file is paired with its record by file name.

**Failed downloads.** beefeater leaves a file that fails to download out of the batch's `downloaded_*.json`, and
logs a warning in that job's `.command.log`. The sample is then missing from the outputs but the run carries on.
If every file in a batch fails, the job fails and is retried by the pipeline's retry strategy.

### 4. Unpacking and publishing (`CRAM_EXTRACT`)

- **Illumina:** `COLLATE_FASTQ` runs `samtools collate | samtools fastq` and writes
  `<ID><raw_reads_suffix>_1.fastq.gz` and `_2.fastq.gz`. Each CRAM is deleted from `work/` once unpacked; a re-run of the pipeline with `-resume` would then have to download that batch from iRODS again, as opposed to reusing cached files.. Set it
  to `false` to keep the CRAMs.
- **ONT fastq:** all `*fastq.gz` files of the collection are concatenated into `<ID>.fastq.gz`.
- **ONT squiggles:** the `.pod5`/`.fast5` files are published to `<outdir>/raw_unbasecalled/<ID>/`.

### Emitted channels

| Channel               | Content                                       |
| --------------------- | --------------------------------------------- |
| `illumina_reads_ch`   | `tuple val(meta), path(read_1), path(read_2)` |
| `ont_reads_ch`        | `tuple val(meta), path(fastq)`                |
| `ont_unbasecalled_ch` | `tuple val(meta), path(pod5/fast5 files)`     |

`meta` holds every iRODS metadata field of the record, including `ID` and `Platform`.

## Avoiding re-downloads

When `prevent_redownloads` is true, `beefeater search` leaves out Illumina records whose reads are already in the
output directory, so they are not downloaded again. A sample counts as present when
`<ID><existing_output_id_suffix><existing_output_extension>` exists:

- `--save_method nested`: under `<outdir>/<ID>/<preexisting_output_dir>/`
- `--save_method flat`: under `<outdir>/<preexisting_output_dir>/`

These four parameters are normally set by the pipeline (in irods_extractor they follow `--only_new_input` and
point at the preprocessed reads when `--preprocessing true`), not by users.

## ONT read selection

For ONT data set `--read_type ont`, and choose with `--ont_get_read_type`:

- `unbasecalled`: raw "squiggles" (`.pod5` for recent data, `.fast5` for older data).
- `reads` (default): the most finalised reads available for each sample, picked by these rules in order:

  1. `offline_basecalled`: `true` ([re-]basecalled with user-requested models) > `false` (standard NPG
     basecalling);
  2. QC status: `pass` > `fail`;
  3. format: `fastq.gz` (basecalled, fully processed) > `bam` (basecalled, needs processing) > `.pod5`
     (unbasecalled, newer format) > `.fast5` (unbasecalled, older format).

  ONT `bam` is not handled by this subworkflow yet (see [Routing by data type](#2-routing-by-data-type-irods_extractor)).

## Parameters

**iRODS query options**

| Option                | Default    | Description                                                                                                 |
| --------------------- | ---------- | ----------------------------------------------------------------------------------------------------------- |
| `--read_type`         | `Illumina` | Platform of the data to retrieve: `Illumina` or `ont`.                                                      |
| `--ont_get_read_type` | `reads`    | (ONT only) `reads` or `unbasecalled`; see [ONT read selection](#ont-read-selection).                        |
| `--search`            | `false`    | Search iRODS and save the metadata only; download nothing.                                                  |
| `--service_user`      | `false`    | Skip beefeater's iRODS login check (`iexpire`). For service users, which authenticate to iRODS another way. |

**iRODS download options**

| Option                        | Default | Description                                                                                                                           |
| ----------------------------- | ------- | ------------------------------------------------------------------------------------------------------------------------------------- |
| `--irods_batch_count`         | `20`    | Maximum number of `BEEFEATER_GET` jobs the non-ONT files are split between. More batches start unpacking sooner but submit more jobs. |
| `--irods_min_batch_size`      | `20`    | Minimum number of non-ONT files in a batch.                                                                                           |
| `--irods_ont_batch_size`      | `1`     | Number of ONT collections per job.                                                                                                    |
| `--irods_max_concurrent_gets` | `3`     | Maximum number of `BEEFEATER_GET` jobs running at once.                                                                               |
| `--irods_pool_size`           | `4`     | iRODS connections held open by each `BEEFEATER_GET` job.                                                                              |
| `--irods_search_pool_size`    | `4`     | iRODS connections held open by the `BEEFEATER_SEARCH` job.                                                                            |
| `--irods_search_tokens`       | `20`    | LSF `seqirods` tokens requested by `BEEFEATER_SEARCH`.                                                                                |
| `--irods_get_tokens`          | `20`    | LSF `seqirods` tokens requested by each `BEEFEATER_GET` job.                                                                          |

**Extraction and output options**

| Option                                         | Default                          | Description                                                                        |
| ---------------------------------------------- | -------------------------------- | ---------------------------------------------------------------------------------- |
| `--save_fastqs`                                | `false` (pipelines may override) | Publish the raw unpacked fastqs.                                                   |
| `--raw_reads_output_dir`                       | `raw_fastqs`                     | Output directory for raw unpacked reads.                                           |
| `--raw_reads_suffix`                           | `_raw`                           | Suffix added to raw fastq names, e.g. `<ID>_raw_1.fastq.gz`.                       |
| `--save_method`                                | `nested`                         | Save outputs in per-sample folders (`nested`) or one folder (`flat`).              |
| `--cleanup_intermediate_files_irods_extractor` | `true`                           | Delete each CRAM from `work/` once unpacked.                                       |
| `--large_data`                                 | `false`                          | Run `COLLATE_FASTQ` from the `week` queue with 10 GB memory, for very large files. |
| `--save_metadata`                              | `true`                           | Save all record metadata to `metadata_irods_queried_<datetime>.csv`.               |
| `--short_metacsv_name`                         | `false`                          | Use a short date instead of a full timestamp in metadata file names.               |

## iRODS connections and LSF tokens

Two separate limits apply to how much data can be queried and downloaded from iRODS (guidance from the Scientific Computing team):

| Resource                   | Maximum Size                    | Controls                                                                               |
| -------------------------- | ------------------------------- | -------------------------------------------------------------------------------------- |
| iRODS database connections | 3000 shared; **10-15 per user** | Live load on the iRODS database. Not enforced: it is up to us.                         |
| LSF `seqirods` tokens      | 1000 per cluster, shared        | Job admission: a job remains pending until the count of requested tokens is available. |

Each connection in a beefeater pool is one iRODS database connection. The search finishes before any download
starts, so the most a run holds at once is `irods_max_concurrent_gets × irods_pool_size`: 12 with the defaults.
`irods_max_concurrent_gets` sets the value of the `maxForks` directive for the `BEEFEATER_GET` process. In general, this directive limits the number of  submitted job for the processes it's applied to. However, this restriction is scoped within one pipeline run only, so two runs by the same user at the same time may submit twice as many jobs, and thus hold up to twice as many connections.
Live connections can be watched at
[the iRODS metrics dashboard](https://metrics.internal.sanger.ac.uk/d/000000187/irods-v4-processes?orgId=1).

## Performance

Study 8279: 1,622 CRAMs, 1.84 TB, download and unpack only (October 2026):

| Version                                                                              | All unpacked    | First sample unpacked |
| ------------------------------------------------------------------------------------ | --------------- | --------------------- |
| irods_extractor v4.0.2 (one beefeater job downloads everything first, 4 connections) | 209.5 min       | 173.5 min             |
| batched, 82 batches of 20 (12 connections)                                           | 76.1 / 78.8 min | 4.7 / 5.5 min         |
| batched by count, 20 batches of ~81 (12 connections)                                 | 73.0 min        | 12.2 min              |

Most of the gain comes from avoiding the bottleneck of waiting for all downloads to complete, instead starting the downstream work (unpacking CRAM files etc.) while later batches are still downloading (asynchronous processing).

## Software versions

| Software              | Version | Image                                             |
| --------------------- | ------- | ------------------------------------------------- |
| Beefeater             | 1.3.0   | `quay.io/sangerpathogens/beefeater:v1.3.0`        |
| Samtools              | 1.17    | `quay.io/biocontainers/samtools:1.17--hd87286a_2` |
| pandas (metadata CSV) | 2.2.1   | `quay.io/sangerpathogens/pandas:2.2.1`            |

Pinned versions are declared in the `modules/*.nf` files.

## Troubleshooting

- **iRODS login error:** run `iinit`. If you see a missing environment file error, check that
  `~/.irods/irods_environment.json` exists, and request an iRODS account if not.
- **No samples found:** check `--studyid`/`--runid`. Use `--search true` to preview the metadata without
  downloading.
- **Samples missing from the outputs:** look in the `.command.log` of the `BEEFEATER_GET` job for that batch for
  `Download for ... failed`.
- **Permission errors:** if a `BEEFEATER_GET` log shows `CAT_NO_ACCESS_PERMISSION`, you do not have access to
  that data. Contact your Scientific Service Representative via dnap-ssr@sanger.ac.uk.
- **"ID(s) are shared by more than one search record":** two iRODS objects would be saved under the same name.
  Narrow the query, or retrieve them with beefeater directly.
- **"iRODS record(s) of a data type this pipeline cannot handle":** those records were left out on purpose; see
  [Routing by data type](#2-routing-by-data-type-irods_extractor).
- **Resuming a failed run:** add `-resume`; completed batches are reused.
