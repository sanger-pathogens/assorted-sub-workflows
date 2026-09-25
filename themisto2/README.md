# themisto2

Nextflow DSL2 sub-workflow library (no `main.nf` of its own) providing `BUILD_COLOUR_INDEX` and `MARKER_FILTERING`, included by parent pipelines such as [lsmd](../../README.md). See that repo's README for the full pipeline documentation -- pipeline steps, parameters, outputs and the `stats.json` field reference.

## Outputs

Every process publishes under `<outdir>/<species>/` ([`modules/publish_paths.nf`](./modules/publish_paths.nf)). Final outputs are always published: `colour_mapping/` (ordered assembly list and label mapping), `index/species_index.thm2`, `index/<group>_marker_index.thm2`, and the candidate-filtering and ATB-check results. Intermediates are published only with `publish_intermediate`: the colour-mapping QC files, and each tool's output under `index/<tool>/<step>/species/` or `.../groups/<group>/`.

## `BUILD_COLOUR_INDEX`

### Inputs

- `samples_ch`: one item per species -- `tuple(meta, metadata, assembly_input)`:
  - `meta.ID` -- species name (output-file prefix / folder name); `MARKER_FILTERING` sets `meta.species` to this value on every per-lineage item it produces, as the join key back to the species-wide outputs
  - `meta.target_groups` -- comma-separated lineage labels to run lineage-specificity filtering for. Empty/absent = every lineage in the metadata with `>= candidate_min_genome_count` genomes (excluding `unclassified`)
  - `metadata` -- that species' metadata table (`.tsv`/`.csv`)
  - `assembly_input` -- a directory of assembly FASTAs, or a `.txt` file listing one assembly path per line (auto-detected)

The parent pipeline builds this channel. In lsmd that's [`subworkflows/manifest_parse.nf`](../../subworkflows/manifest_parse.nf), one item per row of the required `--manifest` TSV (columns `species` / `metadata` / `assemblies` / `target_groups` / `atb_exclude_species`). `--group_label`, `--sample_col` and `--assembly_suffix` stay run-wide params.

### Group-label cleaning

[colour_mapping.py](./bin/colour_mapping.py) reads the metadata as text (so `3` never becomes `3.0`), strips whitespace from headers and the sample and label columns, then applies these rules to each label, in order:

1. **Missing values**: an empty label, or one of `NA`, `N/A`, `#N/A`, `NaN`, `null`, `none`, `unknown`, `missing`, `-`, `?`, `.`, `not applicable`, `not available`, `not collected`, `not provided` (case-insensitive), becomes `unclassified`.
2. **GPSC `;` labels**, only when `--group_label` is `GPSC` (any case). GPS merge-history labels such as `1215;5` or `GPSC1215;GPSC5` become their smallest number, keeping the label's `GPSC` prefix if it has one (`1215;5` → `5`, `GPSC3;28` → `GPSC3`). **Exception:** 235 with 9, in any order or prefix form (`235;9`, `9;235`, `GPSC235;9`, `GPSC9;235`), is kept as its own group `235_9` (`GPSC235_9` with the prefix). It's a mixture of GPSC9 and GPSC235, but current evidence doesn't say to merge the two. A `;` label with a part that isn't a whole number (e.g. `5;abc`) stops the run, listing every bad label. For any other `--group_label`, `;` labels are left as written.

Genomes whose final label is `unclassified` are **always left out of the index**: missing values, labels that already read `unclassified` (any case), and assemblies with no metadata row. Markers are therefore never checked against them. They're listed in `<species>_dropped_unclassified.tsv` (`Sample_ID`, `raw_label`, `reason`: `missing_value` / `labelled_unclassified` / `no_metadata_row`). The run stops if nothing would be left.

`<species>_stats.json` records the missing-value list (`missing_values`), every changed label with its new label, genome count and the rule that changed it (`label_changes`, `changed_by`: `missing_value` or `gpsc_multi`) and `assemblies_dropped_unclassified`.

### Emitted channels

- `sbwt_index`: `tuple(meta, sbwt, lcs)` -- species-wide index.
- `species_export`: `tuple(meta, unitigs, colour_sets, export_metadata, label_mapping)` -- species-wide Themisto2 export; feeds `MARKER_FILTERING`'s `LINEAGE_SPECIFICITY_FILTER` input.
- `checkpoints`: `tuple(meta, row_tsv)` -- per-stage count rows (see "Checkpoint counts" below).

### Lineage-specificity candidate filtering

For each targeted lineage (`meta.target_groups`, or -- when that's blank -- every lineage with `>= candidate_min_genome_count` genomes, excluding `unclassified`), [lineage_specificity_filter.py](./bin/lineage_specificity_filter.py) runs once per species over the **species-wide** export (`export.unitigs.fa` / `export.color_sets.txt` + species `label_mapping.tsv`) and keeps a unitig iff it is:

- lineage-**core** -- present in `>= candidate_min_freq` of the lineage's genomes (`core` 0.95 / `relaxed` 0.5 / `catchall` >0 / literal), and in `>= candidate_min_genome_count` genomes absolutely; and
- lineage-**specific** -- when `specificity_max_outside` is set (default `0.05`), present in `<= specificity_max_outside` of the genomes of *any single other* lineage with `>= candidate_min_genome_count` genomes. `specificity_max_outside = null` disables this and gives the historical core-only behaviour.

Survivors are rebuilt into the candidate index (GGCAT -> SBWT -> Themisto2 build/stats). This replaces the inert `xlin_bg`/`lin_cand` set-diffs (PAT-3570): `sbwt difference` is colour-blind, so a plain set difference can never remove a k-mer a lineage shares with a sister lineage.

## `MARKER_FILTERING`

Everything downstream of the species-wide index build: rebuilds each targeted lineage's candidate markers into their own index, then checks them against ATB for cross-species specificity.

### Inputs

- `species_export_ch`: `tuple(meta, unitigs, colour_sets, export_metadata, label_mapping)` -- `BUILD_COLOUR_INDEX.out.species_export`. `meta.ID` = species id.
- `target_groups_ch`: `tuple(meta, target_groups_string)` -- from the including pipeline's manifest parsing (in lsmd, `MANIFEST_PARSE.out.target_groups`), same slim `[ID: species]` meta as `species_export_ch`, joined in here.
- `atb_target_species_ch`: `tuple(meta, atb_target_species_string)` -- from the including pipeline's manifest parsing (`MANIFEST_PARSE.out.atb_target_species`). In lsmd this is the manifest's `species` value, or a blank string when that name isn't in `--atb_colour_names`, which skips the ATB cross-species check for that species (see below).

### Emitted channels

- `markers`: `tuple(meta, fasta)` -- final candidate markers, `meta.ID` = lineage, `meta.species` set. ATB-checked (`PASS`) for species found in ATB, or the raw rebuilt candidate FASTA (unverified, already logged with a warning) for species that aren't.
- `checkpoints`: `tuple(meta, row_tsv)` -- per-stage count rows (see "Checkpoint counts" below).

### Candidate index rebuild

`LINEAGE_SPECIFICITY_FILTER`'s output (one candidate FASTA per targeted lineage) is wrapped into a colour-list (`CANDIDATE_COLOUR_LIST`) and rebuilt end to end: GGCAT -> SBWT build/check -> Themisto2 build/stats. This is a QC gate only: nothing downstream reads the candidate index. The ATB check below reads GGCAT's candidate unitigs (`GGCAT_CANDIDATE`) directly; the SBWT/Themisto2 rebuild only feeds checkpoints 70/75/80. A lineage whose candidate FASTA comes back empty (nothing cleared the lineage-specificity thresholds) skips the rebuild entirely, with a `log.warn`, rather than failing the run.

> **TODO (backlog):** the ATB check no longer depends on the candidate SBWT/Themisto2 rebuild, so it could go entirely: every count it gives (unitigs, k-mers, strand duplicates) can come from the GGCAT FASTA checkpoint. `THEMISTO2_BUILD` also always publishes the candidate index as `results/<species>/index/<group>_marker_index.thm2`, which reads like a deliverable but is a single-colour index of pre-ATB candidates -- gate it behind `--publish_intermediate` or drop it with the rebuild.

### ATB cross-species check

Replaces the old `bg_excl`/`markers` `sbwt difference` set-diff (PAT-3570: `sbwt difference` is colour-blind and doesn't scale at species-index level). GGCAT's candidate unitigs (`GGCAT_CANDIDATE`, one record per unitig) are pseudoaligned against `ATB-species.thm2` (`THEMISTO2_ATB_PSEUDOALIGN`); [atb_cross_species_filter.py](./bin/atb_cross_species_filter.py) (`FILTER_ATB_MARKERS`) then scores each candidate marker's hit fraction against every ATB species colour and keeps only markers that are solidly within the target species (`>= atb_min_within`, default `0.95`) and essentially absent from every other one (`<= atb_max_outside`, default `0.05`, set independently of `specificity_max_outside`). ATB's `unknown` colour (its catch-all bucket for unassigned/low-confidence genomes) is excluded from the max-outside check entirely by default; add other colours to exclude per species with the manifest's `atb_exclude_species` column (comma-separated; `unknown` is always excluded on top of them).

A species whose manifest `species` name isn't an ATB colour name skips this check (lsmd warns at launch, with the closest ATB names in case it's a typo): its rebuilt candidate markers pass straight through **unchecked**, with a loud `log.warn`, rather than being silently dropped or failing the whole run. Every other species goes through the full check.

Besides the final `PASS` markers, `FILTER_ATB_MARKERS` also writes `FLAG` (off-target leakage, kept for inspection, not forwarded), `ABSENT` (target species never hit at all), a per-marker `validation.tsv`, and a `summary.txt` -- all published under `atb_cross_species/<lineage>/`.

## Dependencies

All software dependencies are containerised (GGCAT, SBWT, Themisto2, and a `pandas` container for [colour_mapping.py](./bin/colour_mapping.py) and [atb_cross_species_filter.py](./bin/atb_cross_species_filter.py)).

## GGCAT `-e` (unitig links)

`GGCAT_CANDIDATE` passes `-e`/`--generate-maximal-unitigs-links` (GGCAT annotates each unitig with its BCALM2-format connectivity links); `GGCAT_SPECIES` does not. This only changes what's written into the FASTA headers (link annotations), not the unitig set GGCAT computes or anything downstream in SBWT/Themisto2 -- confirmed against `ggcat build --help`, which lists it under "Output mode" separately from the actual unitig-generation-mode flags (`--simplitigs`/`--eulertigs`/`--greedy-matchtigs`).

It's candidate-only because it was measured to cost +58% output file size (504MB vs 319MB) on the ~5M-unitig V. cholerae species-wide build, for link data nothing at that scale currently consumes. It's kept on for the candidate rebuild (hundreds-thousands of unitigs, effectively free there) in case a future stitching tool needs it for species whose candidate markers don't self-overlap as cleanly as 7PET's did (PAT-3592).

## Checkpoint counts (`pipeline_counts.tsv`)

Both `BUILD_COLOUR_INDEX` and `MARKER_FILTERING` tap a fixed set of key stages (colour file, GGCAT unitigs, Themisto2 index, exported/dumped FASTA, final markers) through `CHECKPOINT_FASTA`/`CHECKPOINT_THEMISTO` (`modules/checkpoint.nf`, split by input type) as a side channel -- never joined back into the workflow, just counted. Rows from every stage across both subworkflows are combined by the including pipeline's `main.nf` (`collectFile`) into one `pipeline_counts.tsv`, ordered by an `order` key (`BUILD_COLOUR_INDEX` uses 10-40, `MARKER_FILTERING` continues from 50).

Columns, by input `kind`:

| kind | columns populated |
| --- | --- |
| `colourfile` | `n_seqs` (line count) |
| `fasta` | `n_seqs`, `sum_bp`, `min_len`, `median_len`, `max_len` (via `seqkit stats -a`), `n_revcomp_dupes` (via `seqkit rmdup -s`) |
| `themisto` | `n_kmers`, `n_colours`, `n_unitigs` (via `themisto2 stats`) |

**Number of reverse-complement duplicates** (`n_revcomp_dupes`) counts FASTA records that are reverse-complement duplicates of another record already in the same file -- i.e. two records that are the same underlying DNA fragment, just written from opposite strands (`ACGT` vs. its reverse complement `ACGT`->`CGTA`->complemented). A byte-for-byte comparison won't catch these; `seqkit rmdup -s` canonicalises each sequence against its reverse complement before deduping, and compares both strands by default. The rest of this section refers to it by its column name, `n_revcomp_dupes`.

This column exists because `SBWT_DUMP_UNITIGS` (`marker_filtering.nf`, stage `candidate_dumped_fasta`, order 80) reports both strands of every unitig as separate records, so its `n_revcomp_dupes` is expected to be ~100% of `n_seqs` -- this is normal, not a bug. The preceding checkpoint, `candidate_ggcat_unitigs` (order 60, GGCAT's own output before the SBWT round-trip), is expected to show 0 revcomp dupes. Having both stages in `pipeline_counts.tsv` makes that divergence visible on every run without a manual check (found while investigating PAT-3592). The ATB check used to read this dump, which doubled its query and every PASS/FLAG/ABSENT count in `summary.txt`; it now reads GGCAT's output instead.

`candidate_export_unitigs` (order 75, `THEMISTO2_EXPORT_CANDIDATE`) is a checkpoint-only `themisto2 export` of the candidate index, published under `index/themisto2/export/groups/<group>/` and not used by the ATB check. It shows the doubling comes from the candidate SBWT (built with `-r`, so it stores both strands) being read by `sbwt dump-unitigs`, not from Themisto2: the export writes one record per unitig. Checked by hand on the 7PET run of 23 Sep 2026:

| mode | GGCAT (60) | Themisto2 export (75) | SBWT dump (80) | canonical k-mers, all three |
| --- | --- | --- | --- | --- |
| core | 469, 0 dupes | 469, 0 dupes | 938, 469 dupes | 16,106 |
| relaxed | 968, 0 dupes | 968, 0 dupes | 1,936, 968 dupes | 54,447 |
| catchall | 29,206, 0 dupes | 29,206, 0 dupes | 58,410, 29,201 dupes | 2,030,656 |

All three stages hold the identical k-mer set. In catchall, 3 cyclic unitigs (45, 2,119, 2,699 bp) are written from a different start point by each tool, and the dump writes some of them once rather than twice, which is why its dupe count isn't exactly half there. `themisto2 stats`' "Number of forward unitigs (not bidirected)" (checkpoint 70, `n_unitigs`) counts both strands too, so it matches the dump, not the export.

The `seqkit rmdup` dedup check runs unguarded on every `fasta`-kind checkpoint, including species-wide FASTAs (millions of records) -- by design, not oversight. It replaced an older hand-rolled awk canonicalisation that needed a `REVCOMP_CHECK_MAX` size cap because it was too slow to run unguarded at species-wide scale; `seqkit rmdup` doesn't need that cap, benchmarked at ~17s for ~5M records under the fixed `mem_4` label (PAT-3592).
