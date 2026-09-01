# themisto2

Nextflow DSL2 sub-workflow library (no `main.nf` of its own) providing `BUILD_COLOR_INDEX` and `SET_DIFF_CALCULATIONS`, included by parent pipelines such as [lsmd](../../README.md). See that repo's README for the full pipeline documentation -- pipeline steps, parameters, outputs and the `stats.json` field reference.

## `BUILD_COLOR_INDEX`

### Inputs

- `samples_ch`: one item per species -- `tuple(meta, metadata, assembly_input)`:
  - `meta.ID` -- species name (the output folder, and the join key between `species_index` and `lineage_index`)
  - `meta.target_groups` -- comma-separated lineage labels for that species (empty/absent = species-wide only)
  - `metadata` -- that species' metadata table (`.tsv`/`.csv`)
  - `assembly_input` -- a directory of assembly FASTAs, or a `.txt` file listing one assembly path per line

The parent pipeline builds this channel -- e.g. from a run manifest (see lsmd's [`subworkflows/manifest_parse.nf`](../../subworkflows/manifest_parse.nf)). `--group_label`, `--sample_col` and `--assembly_suffix` stay run-wide.

### Emitted channels

`BUILD_COLOR_INDEX`'s public contract is narrower than what it computes internally -- only the SBWT indexes `setdiff_filter.nf` actually needs are exposed:

- `sbwt_index`: `tuple(meta, sbwt, lcs)` -- species-wide index.
- `lineage_index`: `tuple(meta, sbwt, lcs)`, `meta.species` set -- one item per requested lineage (`meta.target_groups`).
- `candidate_index`: `tuple(meta, sbwt, lcs)`, `meta.species` set -- one item per lineage that survived candidate filtering.

The species-wide Themisto2 build/export and `COLOR_MAPPING`'s species `label_mapping`/`stats` still run internally (needed by [lineage_specificity_score.py](./bin/lineage_specificity_score.py), a diagnostic script not yet wired into this subworkflow) -- they're just not part of the emitted contract.

### Dependencies

All software dependencies are containerised (GGCAT, SBWT, Themisto2, and a `pandas` container for [color_mapping.py](./bin/color_mapping.py)).

## `SET_DIFF_CALCULATIONS` (`setdiff_filter.nf`)

Computes step08's index-native set-difference chain (`08_set_difference_filtering.md`, stages `C`/`D`/`F`/`G`). Earlier drafts of this pipeline referred to each stage by that doc's bare letters (`A`-`G`) directly in the code; that's been dropped in favour of descriptive stage/variable names, since bare letters aren't self-documenting to read in `.nf` files without this table open alongside them. This is the current mapping:

| step08 | stage name (code) | formula | produced by | process aliases | emit name |
| --- | --- | --- | --- | --- | --- |
| A | `species_index` | — (built, not diffed) | `BUILD_COLOR_INDEX` (species-wide) | `SBWT_BUILD_SPECIES`, `SBWT_CHECK_SPECIES` | `BUILD_COLOR_INDEX.out.sbwt_index` |
| B | `lineage_index` | — (built, not diffed) | `BUILD_COLOR_INDEX` (per requested lineage) | `SBWT_BUILD_GROUP`, `SBWT_CHECK_GROUP` | `BUILD_COLOR_INDEX.out.lineage_index` |
| C | `bg_excl` | background − `species_index` | `SET_DIFF_CALCULATIONS` | `SBWT_DIFFERENCE_BG_EXCL`, `SBWT_CHECK_BG_EXCL` | `.out.bg_excl` |
| D | `xlin_bg` | `species_index` − `lineage_index` | `SET_DIFF_CALCULATIONS` | `SBWT_DIFFERENCE_XLIN_BG`, `SBWT_CHECK_XLIN_BG` | `.out.xlin_bg` |
| E | `candidate_index` | — (built, not diffed) | `BUILD_COLOR_INDEX` ([core_catchall_filter.py](./bin/core_catchall_filter.py) filtering + GGCAT/SBWT/Themisto2-build-stats rebuild, per lineage) | `SBWT_BUILD_CANDIDATE`, `SBWT_CHECK_CANDIDATE` | `BUILD_COLOR_INDEX.out.candidate_index` |
| F | `lin_cand` | `candidate_index` − `xlin_bg` | `SET_DIFF_CALCULATIONS` | `SBWT_DIFFERENCE_LIN_CAND`, `SBWT_CHECK_LIN_CAND` | `.out.lin_cand` |
| G | `markers` | `lin_cand` − `bg_excl` | `SET_DIFF_CALCULATIONS` (blocked on F) | `SBWT_DIFFERENCE_MARKERS`, `SBWT_CHECK_MARKERS` | `.out.markers` |

`BUILD_COLOR_INDEX.out.lineage_index`/`.out.candidate_index` feed `SET_DIFF_CALCULATIONS`'s `lineage_index_ch`/`candidate_index_ch` in the including pipeline's own `main.nf` (see e.g. [lsmd's main.nf](../../main.nf)).

Every diff (`bg_excl`, `xlin_bg`, `lin_cand`, `markers`) is immediately followed by an `sbwt check` on its output — `sbwt difference` has a confirmed history of silently producing a structurally-corrupt index that LSF still reports as "successfully completed" (see the step08 doc), so never trust a diff's exit code alone.
