# themisto2


## Pipeline Summary

Currently one subworkflow, `BUILD_COLOR_INDEX`:

1. **Colour mapping** ([color_mapping.py](./bin/color_mapping.py)): matches metadata samples to assembly files on disk and writes Themisto's `--file-colors` input, ordered and grouped by `--label_col`. Always writes a species-wide index; optionally also writes one lineage-scoped index per group named in `--target_groups` (Index B -- see [Index B](#index-b-target_groups) below).
2. **GGCAT**: builds unitigs from the species-wide colour file.
3. **SBWT**: builds the SBWT index from the unitigs, then verifies it loads correctly (`sbwt check`, split into its own lighter-weight step).
4. **Themisto2 build**: builds the Themisto2 index from the colour file + verified SBWT index, then sanity-checks it loads (`themisto2 stats`).
5. **Themisto2 export**: exports the index to `export.unitigs.fa`, `export.color_sets.txt` and `export.metadata.txt`.

Steps 2-5 currently only ever run against the species-wide index -- Index B's per-group colour files are written to disk in step 1 but not yet fed through GGCAT/SBWT/Themisto2 (see [Index B](#index-b-target_groups)).

### Parameters

- `--metadata` (required): CSV/TSV file with one row per genome assembly, including a sample-identifier column and a grouping column.
- `--sample_col` (default: `Sample_ID`): metadata column matched against assembly filenames.
- `--label_col` (required): metadata column used to group assemblies into colours (e.g. a GPSC column).
- `--assembly_input` (required): a directory of assembly FASTA files, or a `.txt` file listing one assembly path per line. Which kind it is gets detected automatically.
- `--assembly_suffix` (default: `.contigs.fasta`): suffix appended to `--sample_col` values to form the assembly filename. Only used when `--assembly_input` is a directory.
- `--target_groups` (default: `""`): comma-separated label(s), e.g. `GPSC1,GPSC2`, to also build lineage-scoped indexes for. Leave empty for species-wide only. See [Index B](#index-b-target_groups).
- `--kmer_size` (default: `31`): used consistently across the GGCAT, SBWT and Themisto2 builds. Must match the k-mer size of any existing index this output will later be compared or set-differenced against.
- `--gzip_export` (default: `false`): gzip the Themisto2 export's `color_sets.txt`. Off by default -- compression is single-threaded and slow at full species scale for no benefit beyond disk space.
- `--temp_dir` (default: `""`): shared scratch root for GGCAT/SBWT temp/working dirs. Optional -- falls back to a task-local work directory per tool, which is fine at species scale. Only set this to a Lustre path with confirmed free space when running at full background-DB scale (terabytes free needed). Assumes only one build runs against a given `temp_dir` at a time.

### Inputs

- `metadata_ch`: channel emitting a single metadata file path.
- `assembly_ch`: channel emitting a single assembly directory or paths-list file (matching whichever `--assembly_input` points at).

(Both channels currently only support one metadata/assembly pair per run -- see the `TODO` in [build_color_index.nf](./subworkflows/build_color_index.nf).)

### Outputs

Under `--outdir`:

- `colour_mapping/<ID>/index_species/`: `file_colors_input.txt`, `label_mapping.tsv`, `stats.json` for the species-wide colour set.
- `colour_mapping/<ID>/index_target_group/<group>/`: same three files, scoped to one `--target_groups` label. Only present when `--target_groups` is set.
- `ggcat/<ID>/`: unitigs FASTA built from the species-wide colour file.
- `sbwt/<ID>/`: SBWT index + LCS array built from the unitigs.
- `themisto2/<ID>/build/`: the Themisto2 index (`index.thm2`).
- `themisto2/<ID>/export/`: `export.unitigs.fa`, `export.color_sets.txt` (optionally gzipped), `export.metadata.txt`.

#### stats.json fields

Every `index_species/` and `index_target_group/<group>/` directory gets its own `stats.json`. Some fields only make sense at species-wide scope -- a sample with no label, or an assembly with no metadata row at all, can't be attributed to one specific group, so those fields are simply omitted (not reported as `0`) from per-group `stats.json` files.

| field                                  | species-wide | per-group | meaning                                                                                                                               |
| -------------------------------------- | ------------ | --------- | ------------------------------------------------------------------------------------------------------------------------------------- |
| `index_type`                           | always       | always    | `"species"` or `"target_group"`                                                                                                       |
| `target_group`                         | always       | always    | `"species_wide"`, or the group's label                                                                                                |
| `metadata_column`                      | always       | always    | the `--label_col` value used                                                                                                          |
| `samples_dropped_missing_label`        | always       | omitted   | samples with no value in `--label_col` at all                                                                                         |
| `samples_dropped_missing_assembly`     | always       | always    | samples with a label but no matching assembly on disk -- this one genuinely varies by group, so per-group indexes get their own count |
| `assemblies_excluded_missing_metadata` | always       | omitted   | assembly files on disk with no matching metadata row                                                                                  |
| `total_assemblies_written`             | always       | always    | assemblies actually written to this index's `file_colors_input.txt`                                                                   |

### Index B (`--target_groups`)

Index A (species-wide) scales poorly once a species has many groups -- s_pneumoniae alone has 765+ GPSCs. `--target_groups` lets you additionally build colour files scoped to just the lineages you care about, without re-running the whole pipeline per lineage.

As it stands today, setting `--target_groups` does not skip or shrink the species-wide build -- `index_species/` is always built regardless. The per-group colour files under `index_target_group/<group>/` are written by [color_mapping.py](./bin/color_mapping.py) and emitted by `COLOR_MAPPING` (`target_group_indexes`, optional), but nothing downstream consumes them yet -- GGCAT, SBWT and Themisto2 still only ever build from `index_species/`. Wiring a per-group build/export path is separate, not-yet-done work.

### Dependencies

All software dependencies are containerised (GGCAT, SBWT, Themisto2, and a `pandas` container for [color_mapping.py](./bin/color_mapping.py)).

## Set-difference subworkflow (`setdiff_filter.nf`)

Computes step08's index-native set-difference chain
(`08_set_difference_filtering.md`, stages `C`/`D`/`F`/`G`) as workflow
`SET_DIFF_CALCULATIONS`. Earlier drafts of this pipeline referred to each
stage by that doc's bare letters (`A`-`G`) directly in the code; that's been
dropped in favour of descriptive stage/variable names, since bare letters
aren't self-documenting to read in `.nf` files without this table open
alongside them. This is the current mapping:

| step08 | stage name (code) | formula | produced by | process aliases | emit name |
| --- | --- | --- | --- | --- | --- |
| A | `species_index` | — (built, not diffed) | `BUILD_COLOR_INDEX` (species-wide `index_species/`) | `SBWT_BUILD`, `SBWT_CHECK` | `BUILD_COLOR_INDEX.out.sbwt_index` |
| B | `lineage_index` | — (built, not diffed) | **not built yet** — `--target_groups` colour files already exist (`index_target_group/<group>/`, see [Index B](#index-b-target_groups) above) but nothing downstream builds an SBWT index from them yet | n/a | n/a |
| C | `bg_excl` | background − `species_index` | `SET_DIFF_CALCULATIONS` | `SBWT_DIFFERENCE_BG_EXCL`, `SBWT_CHECK_BG_EXCL` | `.out.bg_excl` |
| D | `xlin_bg` | `species_index` − `lineage_index` | `SET_DIFF_CALCULATIONS` (blocked on B) | `SBWT_DIFFERENCE_XLIN_BG`, `SBWT_CHECK_XLIN_BG` | `.out.xlin_bg` |
| E | `candidate_index` | — (built, not diffed) | **not implemented anywhere yet** — Step07-derived per-lineage candidate core | n/a | n/a |
| F | `lin_cand` | `candidate_index` − `xlin_bg` | `SET_DIFF_CALCULATIONS` (blocked on E and D) | `SBWT_DIFFERENCE_LIN_CAND`, `SBWT_CHECK_LIN_CAND` | `.out.lin_cand` |
| G | `markers` | `lin_cand` − `bg_excl` | `SET_DIFF_CALCULATIONS` (blocked on F) | `SBWT_DIFFERENCE_MARKERS`, `SBWT_CHECK_MARKERS` | `.out.markers` |

Every diff (`bg_excl`, `xlin_bg`, `lin_cand`, `markers`) is immediately
followed by an `sbwt check` on its output — `sbwt difference` has a
confirmed history of silently producing a structurally-corrupt index that
LSF still reports as "successfully completed" (see the step08 doc), so
never trust a diff's exit code alone.
