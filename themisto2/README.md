# themisto2

Nextflow DSL2 sub-workflow library (no `main.nf` of its own) providing `BUILD_COLOR_INDEX` and `SET_DIFF_CALCULATIONS`, included by parent pipelines such as [lsmd](../../README.md). See that repo's README for the full pipeline documentation -- pipeline steps, parameters, outputs and the `stats.json` field reference.

## `BUILD_COLOR_INDEX`

### Inputs

- `samples_ch`: one item per species -- `tuple(meta, metadata, assembly_input)`:
  - `meta.ID` -- species name (the output folder, and the join key between `species_index` and `candidate_index`)
  - `meta.target_groups` -- comma-separated lineage labels to run step07 lineage-specificity filtering for (empty/absent = species-wide index only)
  - `metadata` -- that species' metadata table (`.tsv`/`.csv`)
  - `assembly_input` -- a directory of assembly FASTAs, or a `.txt` file listing one assembly path per line

The parent pipeline builds this channel -- e.g. from a run manifest (see lsmd's [`subworkflows/manifest_parse.nf`](../../subworkflows/manifest_parse.nf)). `--group_label`, `--sample_col` and `--assembly_suffix` stay run-wide.

### Emitted channels

`BUILD_COLOR_INDEX`'s public contract is narrower than what it computes internally -- only the SBWT indexes `setdiff_filter.nf` actually needs are exposed:

- `sbwt_index`: `tuple(meta, sbwt, lcs)` -- species-wide index.
- `candidate_index`: `tuple(meta, sbwt, lcs)`, `meta.species` set -- one item per lineage in `meta.target_groups` that survived step07 lineage-specificity filtering.

The species-wide Themisto2 build/export and `COLOR_MAPPING`'s species `label_mapping`/`stats` still run internally -- step07 ([lineage_specificity_filter.py](./bin/lineage_specificity_filter.py)) reads them -- they're just not part of the emitted contract.

### Step07 -- lineage-specificity candidate filtering

For each lineage in `meta.target_groups`, [lineage_specificity_filter.py](./bin/lineage_specificity_filter.py) runs once per species over the **species-wide** export (`export.unitigs.fa` / `export.color_sets.txt` + species `label_mapping.tsv`) and keeps a unitig iff it is:

- lineage-**core** -- present in `>= candidate_min_freq` of the lineage's genomes (`core` 0.95 / `relaxed` 0.5 / `catchall` >0 / literal), and in `>= candidate_min_genome_count` genomes absolutely; and
- lineage-**specific** -- when `specificity_max_outside` is set (default `0.05`), present in `<= specificity_max_outside` of the genomes of *any single other* lineage with `>= candidate_min_genome_count` genomes. `specificity_max_outside = null` disables this and gives the historical core-only behaviour.

Survivors are rebuilt into the candidate index `E` (GGCAT -> SBWT -> Themisto2 build/stats). This replaces the inert `xlin_bg`/`lin_cand` set-diffs (PAT-3570): `sbwt difference` is colour-blind, so a plain set difference can never remove a k-mer a lineage shares with a sister lineage.

### Dependencies

All software dependencies are containerised (GGCAT, SBWT, Themisto2, and a `pandas` container for [color_mapping.py](./bin/color_mapping.py)).

## `SET_DIFF_CALCULATIONS` (`setdiff_filter.nf`)

Computes step08's index-native set differences (`08_set_difference_filtering.md`). Earlier drafts referred to each stage by that doc's bare letters (`A`-`G`) directly in the code; that's been dropped in favour of descriptive stage/variable names. The `xlin_bg` (`D`) and `lin_cand` (`F`) cross-lineage set-diffs were removed in PAT-3570 -- `sbwt difference` is colour-blind, so `E − (A − B) == E` whenever `E`'s k-mers are lineage-core, which they always are. Cross-lineage specificity is now the differential-frequency filter at step07. Current mapping:

| step08 | stage name (code) | formula | produced by | process aliases | emit name |
| --- | --- | --- | --- | --- | --- |
| A | `species_index` | — (built, not diffed) | `BUILD_COLOR_INDEX` (species-wide) | `SBWT_BUILD_SPECIES`, `SBWT_CHECK_SPECIES` | `BUILD_COLOR_INDEX.out.sbwt_index` |
| C | `bg_excl` | background − `species_index` | `SET_DIFF_CALCULATIONS` | `SBWT_DIFFERENCE_BG_EXCL`, `SBWT_CHECK_BG_EXCL` | `.out.bg_excl` |
| E | `candidate_index` | — (built, not diffed) | `BUILD_COLOR_INDEX` (step07 [lineage_specificity_filter.py](./bin/lineage_specificity_filter.py) + GGCAT/SBWT/Themisto2-build-stats rebuild, per lineage) | `SBWT_BUILD_CANDIDATE`, `SBWT_CHECK_CANDIDATE` | `BUILD_COLOR_INDEX.out.candidate_index` |
| G | `markers` | `candidate_index` − `bg_excl` | `SET_DIFF_CALCULATIONS` | `SBWT_DIFFERENCE_MARKERS`, `SBWT_CHECK_MARKERS` | `.out.markers` |

`BUILD_COLOR_INDEX.out.sbwt_index`/`.out.candidate_index` feed `SET_DIFF_CALCULATIONS`'s `species_index_ch`/`candidate_index_ch` in the including pipeline's own `main.nf`.

Every diff (`bg_excl`, `markers`) is immediately followed by an `sbwt check` on its output — `sbwt difference` has a confirmed history of silently producing a structurally-corrupt index that LSF still reports as "successfully completed" (see the step08 doc), so never trust a diff's exit code alone.
