# themisto2


## Pipeline Summary

Currently one subworkflow, `BUILD_COLOR_INDEX`:

1. **Colour mapping** ([color_mapping.py](./bin/color_mapping.py)): matches metadata samples to assembly files on disk and writes Themisto's `--file-colors` input, ordered and grouped by `--group_label`. Always writes a species-wide index (Index A); optionally also writes one lineage-scoped index per group named in `--target_groups` (Index B -- see [Index B](#index-b-target_groups) below).
2. **GGCAT**: builds unitigs from the colour file.
3. **SBWT**: builds the SBWT index from the unitigs, then verifies it loads correctly (`sbwt check`, split into its own lighter-weight step).
4. **Themisto2 build**: builds the Themisto2 index from the colour file + verified SBWT index, then sanity-checks it loads (`themisto2 stats`).
5. **Themisto2 export**: exports the index to `export.unitigs.fa`, `export.color_sets.txt` and `export.metadata.txt`.

Steps 2-5 run for both Index A (species-wide) and Index B (per-lineage, when `--target_groups` is set) as two separate parallel pipelines -- same steps, same processes, just aliased (`*_SPECIES` / `*_GROUP`) so each can be invoked once per Nextflow workflow scope.

6. **Candidate index filtering** ([core_catchall_filter.py](./bin/core_catchall_filter.py), per lineage, only when `--target_groups` is set): filters that lineage's own Step 5 export down to a candidate marker unitig FASTA using `--candidate_min_freq`/`--candidate_min_genome_count`. This is Python-derived, not yet a real index.
7. **Candidate index rebuild**: the filtered FASTA is rebuilt through GGCAT → SBWT build/check (a real index, Index E) → Themisto2 build/stats (QC gate only -- confirms the rebuild is structurally sound, no export, since nothing downstream reads it). If nothing survives filtering, the FASTA is empty and this rebuild is skipped for that lineage (logged via `log.warn`).

### Parameters

- `--metadata` (required): CSV/TSV file with one row per genome assembly, including a sample-identifier column and a grouping column.
- `--sample_col` (default: `Sample_ID`): metadata column matched against assembly filenames.
- `--group_label` (required): metadata column used to group assemblies into colours (e.g. a GPSC column).
- `--assembly_input` (required): a directory of assembly FASTA files, or a `.txt` file listing one assembly path per line. Which kind it is gets detected automatically.
- `--assembly_suffix` (default: `.contigs.fasta`): suffix appended to `--sample_col` values to form the assembly filename. Only used when `--assembly_input` is a directory.
- `--target_groups` (default: `""`): comma-separated label(s), e.g. `GPSC1,GPSC2`, to also build lineage-scoped indexes for. Leave empty for species-wide only. See [Index B](#index-b-target_groups).
- `--kmer_size` (default: `31`): used consistently across the GGCAT, SBWT and Themisto2 builds. Must match the k-mer size of any existing index this output will later be compared or set-differenced against.
- `--gzip_export` (default: `false`): gzip the Themisto2 export's `color_sets.txt`. Off by default -- compression is single-threaded and slow at full species scale for no benefit beyond disk space.
- `--temp_dir` (default: `""`): scratch root for GGCAT/SBWT temp/working dirs, per meta.ID subdirectory (`<temp_dir>/<tool>/<meta.ID>/`, so concurrent species/lineage/candidate tasks don't collide). Optional -- falls back to a task-local work directory per tool, which is fine at species scale. Only set this to a Lustre path with confirmed free space when running at full background-DB scale (terabytes free needed).
- `--candidate_min_freq` (default: `"core"`): presence-fraction threshold for candidate index filtering (Step 6), only used when `--target_groups` is set. Pass a named preset (`core` >=0.95, `relaxed` >=0.5, `catchall` present in >=1 genome) or your own fraction directly, e.g. `"0.8"`. See [core_catchall_filter.py](./bin/core_catchall_filter.py)'s own docstring for the exact semantics.
- `--candidate_min_genome_count` (default: `5`): absolute genome-count floor for candidate index filtering, applied alongside `--candidate_min_freq` (a unitig needs both to be kept).

### Inputs

- `metadata_ch`: channel emitting a single metadata file path.
- `assembly_ch`: channel emitting a single assembly directory or paths-list file (matching whichever `--assembly_input` points at).

(Both channels currently only support one metadata/assembly pair per run -- see the `TODO` in [build_color_index.nf](./subworkflows/build_color_index.nf).)

### Outputs

Under `--outdir`, `<ID>` is either the species run's own ID (Index A) or a `--target_groups` label like `GPSC1` (Index B and its candidate rebuild) -- both share the same directory shape:

- `colour_mapping/<ID>/index_species/`: `species_file_colors_input.txt`, `species_label_mapping.tsv`, `species_stats.json` for the species-wide colour set.
- `colour_mapping/<ID>/index_target_group/<group>/`: same three files, tagged `<group>_*` instead of `species_*` and scoped to one `--target_groups` label. Only present when `--target_groups` is set.
- `ggcat/<ID>/`: unitigs FASTA built from the colour file.
- `sbwt/<ID>/`: SBWT index + LCS array built from the unitigs.
- `themisto2/<ID>/build/`: the Themisto2 index (`index.thm2`).
- `themisto2/<ID>/export/`: `export.unitigs.fa`, `export.color_sets.txt` (optionally gzipped), `export.metadata.txt`. Not written for the candidate rebuild (Step 7) -- that stage only builds+stats as a QC gate, never exports.
- `candidate_filter/<lineage>/`: `{lineage_id}_{min_freq_label}_candidate_unitigs.fasta` + `_stats.txt` from Step 6, before the Step 7 rebuild. Only present when `--target_groups` is set.

#### Emitted channels

`BUILD_COLOR_INDEX`'s public contract is narrower than what it computes internally -- only the SBWT indexes `set_diff_calculations.nf` actually needs are exposed:

- `sbwt_index`: `tuple(meta, sbwt, lcs)` -- Index A (species-wide).
- `lineage_index`: `tuple(meta, sbwt, lcs)`, `meta.species` set -- Index B, one item per `--target_groups` lineage.
- `candidate_index`: `tuple(meta, sbwt, lcs)`, `meta.species` set -- Index E (the Step 7 rebuild), one item per lineage that survived Step 6's filtering.

The species-wide Themisto2 build/export and `COLOR_MAPPING`'s species `label_mapping`/`stats` still run internally (needed by [lineage_specificity_score.py](./bin/lineage_specificity_score.py), a diagnostic script not yet wired into this subworkflow) -- they're just not part of the emitted contract.

#### stats.json fields

Every `index_species/` and `index_target_group/<group>/` directory gets its own `stats.json`. Some fields only make sense at species-wide scope -- a sample with no label, or an assembly with no metadata row at all, can't be attributed to one specific group, so those fields are simply omitted (not reported as `0`) from per-group `stats.json` files.

| field                                  | species-wide | per-group | meaning                                                                                                                               |
| -------------------------------------- | ------------ | --------- | ------------------------------------------------------------------------------------------------------------------------------------- |
| `index_type`                           | always       | always    | `"species"` or `"target_group"`                                                                                                       |
| `target_group`                         | always       | always    | `"species_wide"`, or the group's label                                                                                                |
| `metadata_column`                      | always       | always    | the `--group_label` value used                                                                                                        |
| `samples_dropped_missing_label`        | always       | omitted   | samples with no value in `--group_label` at all                                                                                       |
| `samples_dropped_missing_assembly`     | always       | always    | samples with a label but no matching assembly on disk -- this one genuinely varies by group, so per-group indexes get their own count |
| `assemblies_excluded_missing_metadata` | always       | omitted   | assembly files on disk with no matching metadata row                                                                                  |
| `total_assemblies_written`             | always       | always    | assemblies actually written to this index's `file_colors_input.txt`                                                                   |

### Index B (`--target_groups`)

Index A (species-wide) scales poorly once a species has many groups -- s_pneumoniae alone has 765+ GPSCs. `--target_groups` lets you additionally build colour files scoped to just the lineages you care about, without re-running the whole pipeline per lineage.

Setting `--target_groups` does not skip or shrink the species-wide build -- `index_species/` is always built regardless. The per-group colour files under `index_target_group/<group>/` are written by [color_mapping.py](./bin/color_mapping.py) and fanned out (one Nextflow channel item per lineage, tagged `meta = [ID: <group>, species: <species run's ID>]`) into the same GGCAT → SBWT build/check → Themisto2 build/stats/export chain Index A uses, via aliased `*_GROUP` processes. Each lineage's own export then feeds Step 6/7's candidate filtering + rebuild.

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
| A | `species_index` | — (built, not diffed) | `BUILD_COLOR_INDEX` (species-wide) | `SBWT_BUILD_SPECIES`, `SBWT_CHECK_SPECIES` | `BUILD_COLOR_INDEX.out.sbwt_index` |
| B | `lineage_index` | — (built, not diffed) | `BUILD_COLOR_INDEX` (per `--target_groups` lineage, see [Index B](#index-b-target_groups) above) | `SBWT_BUILD_GROUP`, `SBWT_CHECK_GROUP` | `BUILD_COLOR_INDEX.out.lineage_index` |
| C | `bg_excl` | background − `species_index` | `SET_DIFF_CALCULATIONS` | `SBWT_DIFFERENCE_BG_EXCL`, `SBWT_CHECK_BG_EXCL` | `.out.bg_excl` |
| D | `xlin_bg` | `species_index` − `lineage_index` | `SET_DIFF_CALCULATIONS` | `SBWT_DIFFERENCE_XLIN_BG`, `SBWT_CHECK_XLIN_BG` | `.out.xlin_bg` |
| E | `candidate_index` | — (built, not diffed) | `BUILD_COLOR_INDEX` ([core_catchall_filter.py](./bin/core_catchall_filter.py) filtering + GGCAT/SBWT/Themisto2-build-stats rebuild, per lineage) | `SBWT_BUILD_CANDIDATE`, `SBWT_CHECK_CANDIDATE` | `BUILD_COLOR_INDEX.out.candidate_index` |
| F | `lin_cand` | `candidate_index` − `xlin_bg` | `SET_DIFF_CALCULATIONS` | `SBWT_DIFFERENCE_LIN_CAND`, `SBWT_CHECK_LIN_CAND` | `.out.lin_cand` |
| G | `markers` | `lin_cand` − `bg_excl` | `SET_DIFF_CALCULATIONS` (blocked on F) | `SBWT_DIFFERENCE_MARKERS`, `SBWT_CHECK_MARKERS` | `.out.markers` |

`SET_DIFF_CALCULATIONS` itself isn't called from anywhere yet -- no `main.nf` exists in this repo (a sub-workflow library, included by a parent pipeline), so B and E now having real producers means the *data* is available, but wiring `BUILD_COLOR_INDEX.out.lineage_index`/`.out.candidate_index` into `SET_DIFF_CALCULATIONS`'s `lineage_index_ch`/`candidate_index_ch` is still a `main.nf`-level task, not yet done. See the TODO block at the bottom of [setdiff_filter.nf](./subworkflows/setdiff_filter.nf).

Every diff (`bg_excl`, `xlin_bg`, `lin_cand`, `markers`) is immediately
followed by an `sbwt check` on its output — `sbwt difference` has a
confirmed history of silently producing a structurally-corrupt index that
LSF still reports as "successfully completed" (see the step08 doc), so
never trust a diff's exit code alone.
