# Strain Mapper subworkflow

Maps paired short reads to a reference genome, calls variants and generates a consensus sequence. Used by the [strain_mapper](https://gitlab.internal.sanger.ac.uk/sanger-pathogens/pipelines/strain_mapper) pipeline.

## Multiple reference support

The subworkflow enables mapping every sample in a run against one reference genome, supplied with `--reference`. It also accepts a **reference manifest**, which assigns a reference to each sample individually. A single run can therefore mix references — useful when the input spans several species or strains, or when each isolate should be mapped against its own closest assembly.

### Supplying references

| Option                 | Description                                                                 |
| ---------------------- | --------------------------------------------------------------------------- |
| `--reference`          | A single reference FASTA, used for any sample without a manifest entry.      |
| `--reference_manifest` | A CSV assigning a reference FASTA per sample ID.                             |

At least one of the two must be supplied; the pipeline prints the help message and exits if neither is given. They are intended to be used together: the manifest covers the samples that need a specific reference, and `--reference` supplies the reference genome for the rest of samples without specific reference allocation.

### Reference manifest format

A CSV with the required header `ID,reference`:

```
ID,reference
sampleA,/path/to/strain_1.fasta
sampleB,/path/to/strain_2.fasta
sampleC,NA
```

- **`ID`** must match a sample ID from the reads input — the `ID` column of `--manifest_of_reads`, or the ID that `mixed_input` derives for iRODS, ENA and directory input.
- **`reference`** is a path to a FASTA file. Paths are validated while the manifest is parsed, and the run fails immediately with the offending path if a file is missing.
- **`NA`** means "no specific reference for this sample". The row is dropped and the sample falls back to `--reference` (note that it is not required to document samples that do not have a specific reference, they can just be omitted from the manifest of references).

### How each sample gets its reference

For every sample emerging from `mixed_input`:

1. If its ID appears in the reference manifest with a real path, that reference is used.
2. Otherwise `--reference` is used.
3. If neither applies:
  a. if `--drop_without_ref` is `true`: the reference-less sample is **dropped from the run**
  b. if `--drop_without_ref` is `false` (default): an error will be raised, interrupting the pipeline.

If `--drop_without_ref` is `true`, running with `--reference_manifest` alone and a manifest that does not cover every sample may therefore produce fewer results than samples submitted (a warning is issued). Supply `--reference` as a fallback unless you deliberately want to restrict the run to the manifested samples.

### Indexing

Indexes are built **once per distinct reference path**, not once per sample, so a reference shared by many samples is indexed a single time.

For each reference, an existing index sitting next to the FASTA is reused if it is complete; otherwise the subworkflow builds one. "Complete" is strict — a partial set is ignored and rebuilt rather than passed to the aligner:

| Index                 | Reused when these files exist next to the reference                                       |
| --------------------- | ----------------------------------------------------------------------------------------- |
| Bowtie2 (`--mapper bowtie2`) | all six of `<basename>.{1,2,3,4,rev.1,rev.2}.bt2`, where `<basename>` is the FASTA name without its extension (e.g. `strain_1.1.bt2` for `strain_1.fasta`) |
| BWA (`--mapper bwa`)  | all five of `<reference>.{amb,ann,bwt,pac,sa}` (e.g. `strain_1.fasta.amb`)                 |
| Samtools faidx        | `<reference>.fai`                                                                           |

Note the Bowtie2 naming: the index prefix is the FASTA name **minus its extension**, matching what this subworkflow writes when it builds one. An index built elsewhere as `bowtie2-build strain_1.fasta strain_1.fasta` (producing `strain_1.fasta.1.bt2`) is not recognised and will be rebuilt.

### Output

Per-sample output is unchanged and still lands under `${outdir}/<sample_ID>/`. The consensus FASTA carries the reference it was called against in its filename:

```
results/<sample_ID>/curated_consensus/<sample_ID>_<reference_basename>.fa
```

Indexes built by the subworkflow are published flat, to `${outdir}/bowtie2`, `${outdir}/bwa` and `${outdir}/sorted_ref`.

## Workflow interface

The `STRAIN_MAPPER` workflow takes a single channel of per-sample tuples, each already carrying its own reference:

```
[ meta, reads_1, reads_2, reference ]
```

This replaces the previous two-input form, `STRAIN_MAPPER(ch_reads, reference)`, where the reference was a separate value channel broadcast to every sample. Reference resolution — joining the reads to the reference manifest and applying the `--reference` fallback — happens in the calling pipeline's `main.nf`, before `STRAIN_MAPPER` is invoked.

Internally, `INDEX_REF` emits three channels, each keyed on the original reference path as a string so that indexes can be joined back to the reads that need them:

```
ch_bt2_index    // [ ref_key, reference, bt2_index_files ]
ch_bwa_index    // [ ref_key, reference, bwa_index_files ]
ch_ref_index    // [ ref_key, reference, faidx ]
```

The string key is necessary because an index process re-emits its reference as a work directory path, which no longer compares equal to the path the reads carry.

`REF_MANIFEST_PARSE` is a separate entry point that parses and validates the reference manifest, emitting `[ meta, reference ]` with `NA` rows already removed.

## Parameters

```
 Read mapping input parameters
      --reference
            default: null
            Path to reference genome in fasta format. At least one of this or
            --reference_manifest needs to be supplied.
      --reference_manifest
            default: null
            Path to a manifest listing reference genomes in fasta format. Required
            headers are ID,reference. Input samples not listed in the reference
            manifest are mapped against the reference provided with --reference.
```

See `--help` for the full parameter list, including mapping and variant filtering options.

## Known limitations

- **Index publish directories are flat.** Two references with the same filename in different directories overwrite each other's published index files under `${outdir}/bowtie2`, `${outdir}/bwa` and `${outdir}/sorted_ref`. Mapping itself is unaffected, as indexes are keyed on the full reference path.
- **Reference manifest IDs are not cross-checked against the reads input.** A row whose `ID` matches no sample is not reported. This is intentional as it enables maintaining a large sample-to-reference manifest or "database" that spans many species and strains, not all of which are meant to be included in a single pipeline run.
- **Samples with no reference are dropped, with a warning at the beginning of the workflow**.
