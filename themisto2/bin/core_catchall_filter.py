#!/usr/bin/env python3

"""
Filters a lineage-scoped Themisto2 export (one lineage's own already-built
index -- see color_mapping.py's --target_groups) down to a candidate
marker unitig set, written as a FASTA for GGCAT to rebuild.

A unitig's within-lineage presence fraction is size/n_colors, read
straight off export.color_sets.txt's 'size=M' field -- no species-wide
comparison needed, since every colour here already is a lineage genome.
n_colors comes from this export's own export.metadata.txt ('num_colors=');
the lineage ID comes from color_mapping.py's label_mapping.tsv.

A unitig is kept only if it clears BOTH:
  - --min-freq: fraction >= min_freq (fraction > 0 if min_freq is 0.0, so
    nothing totally absent from the lineage ever passes). Pass a preset
    ('core'=0.95, 'relaxed'=0.5, 'catchall'=0.0) or your own value.
  - --min-genome-count (default 5): an absolute genome-count floor, so a
    thin fraction isn't satisfied by too few genomes to trust.

Output: '{lineage_id}_{min_freq_label}_candidate_unitigs.fasta', empty
(with a warning) if nothing passes -- BUILD_COLOR_INDEX then skips the
rebuild for that lineage.
"""

import argparse
import gzip
import sys
from pathlib import Path
import numpy as np

##############################################################################
### Presets
#
# Named shortcuts for the one real parameter, the presence-fraction cutoff
# (fraction >= min_freq). --min-freq accepts either one of these names or
# a custom fraction written directly.

PRESET_MIN_FREQ = {
    "core": 0.95,
    "relaxed": 0.5,
    "catchall": 0.0,
}

DEFAULT_MIN_FREQ_PRESET = "core"
DEFAULT_MIN_GENOME_COUNT = 5

##############################################################################
### CLI

def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute within-lineage presence for unitigs in a "
                     "lineage-scoped Themisto2 export (BUILD_COLOR_INDEX's "
                     "--target_groups output), then apply --min-freq/"
                     "--min-genome-count to write a filtered candidate "
                     "marker set as a FASTA, ready for GGCAT to rebuild.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--unitigs", required=True, type=Path,
                         help="export.unitigs.fa(.gz) from this lineage's own "
                              "Themisto2 export (BUILD_COLOR_INDEX's lineage "
                              "branch) -- NOT the species-wide export.")
    parser.add_argument("--color-sets", required=True, type=Path,
                         help="export.color_sets.txt(.gz) from the same "
                              "lineage-scoped Themisto2 export. Format per "
                              "line: 'color_set_id=N size=M c1 c2 c3 ...' -- "
                              "only id=N and size=M are read.")
    parser.add_argument("--export-metadata", required=True, type=Path,
                         help="This lineage's own export.metadata.txt from "
                              "the same Themisto2 export ('key=value' "
                              "lines) -- supplies n_colors from its "
                              "num_colors= line, Themisto2's own "
                              "authoritative genome count for this index.")
    parser.add_argument("--label-mapping", required=True, type=Path,
                         help="This lineage's own "
                              "index_target_group/<group>/label_mapping.tsv, "
                              "already written by color_mapping.py -- "
                              "supplies the lineage ID from its 'label' "
                              "column (one value, shared by every row).")
    parser.add_argument("--output-dir", required=True, type=Path,
                         help="Directory to write the filtered candidate "
                              "unitigs FASTA to, named "
                              "'{lineage_id}_{min_freq_label}_candidate_unitigs.fasta'. "
                              "Empty (no records) if nothing clears both "
                              "thresholds. Created if it doesn't exist.")
    parser.add_argument("--stats-output-dir", type=Path, default=None,
                         help="Optional directory to write a summary stats "
                              "file (total/kept/dropped counts), named "
                              "'{lineage_id}_{min_freq_label}_stats.txt'. "
                              "Created if it doesn't exist. If omitted, no "
                              "stats are written.")
    parser.add_argument("-f", "--min-freq",
                         default=DEFAULT_MIN_FREQ_PRESET,
                         help="Minimum fraction of the lineage's genomes a "
                              "unitig must be present in to be kept, e.g. "
                              "0.95 = present in at least 95%% of genomes. "
                              f"Named presets: 'core' "
                              f"({PRESET_MIN_FREQ['core']}), 'relaxed' "
                              f"({PRESET_MIN_FREQ['relaxed']}), 'catchall' "
                              "(present in at least 1 genome) -- or set "
                              "your own number directly, e.g. --min-freq "
                              f"0.8. Default: '{DEFAULT_MIN_FREQ_PRESET}'.")
    parser.add_argument("--min-genome-count",
                         type=int,
                         default=DEFAULT_MIN_GENOME_COUNT,
                         help="Minimum absolute number of genomes a unitig "
                              "must be present in to be kept, applied ON "
                              "TOP OF --min-freq (a unitig needs both). "
                              "Guards against the fraction alone being "
                              "satisfied by too few genomes -- e.g. a "
                              "unitig found in just 1 genome still "
                              f"clearing min_freq=0.0. Default "
                              f"{DEFAULT_MIN_GENOME_COUNT}.")
    return parser.parse_args()

##############################################################################
### I/O helpers

def _open(path):
    path = Path(path)
    if path.exists():
        return gzip.open(path, "rt") if path.suffix == ".gz" else open(path)
    gz_path = Path(str(path) + ".gz")
    if gz_path.exists():
        return gzip.open(gz_path, "rt")
    no_gz_path = Path(str(path)[:-3]) if str(path).endswith(".gz") else None
    if no_gz_path and no_gz_path.exists():
        return open(no_gz_path)
    raise FileNotFoundError(f"Could not find file: {path} (also tried .gz variant)")


def _open_out(path):
    path = Path(path)
    return gzip.open(path, "wt") if path.suffix == ".gz" else open(path, "w")

##############################################################################
### I/O loading + streaming

def load_n_colors(export_metadata_path) -> int:
    # Themisto2's own authoritative genome count for this index -- not recomputed here.
    with _open(export_metadata_path) as fh:
        for line in fh:
            key, _, value = line.strip().partition("=")
            if key == "num_colors":
                return int(value)
    raise ValueError(f"{export_metadata_path} has no 'num_colors=' line")


def load_lineage_id(label_mapping_path) -> str:
    # Every row shares the same label (file is already scoped to one target_group).
    with _open(label_mapping_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        label_col = header.index("label")
        first_row = fh.readline()
        if not first_row:
            raise ValueError(f"{label_mapping_path} has no data rows -- can't determine lineage_id")
        return first_row.rstrip("\n").split("\t")[label_col]


def load_colorset_sizes(color_sets_path) -> np.ndarray:
    # Dense array of size=M, indexed by color_set_id. c1 c2 c3 ... is never parsed --
    # every colour in a lineage-scoped export already IS a lineage genome, so size=M
    # already IS the within-lineage genome count.
    colorset_ids = []
    sizes = []
    with _open(color_sets_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            id_field, size_field = line.split(maxsplit=2)[:2]
            colorset_ids.append(int(id_field.split("=")[1]))
            sizes.append(int(size_field.split("=")[1]))

    if not colorset_ids:
        return np.zeros(0, dtype=np.int64)

    colorset_ids_arr = np.array(colorset_ids, dtype=np.int64)
    sizes_arr = np.array(sizes, dtype=np.int64)
    dense_sizes = np.zeros(int(colorset_ids_arr.max()) + 1, dtype=np.int64)
    dense_sizes[colorset_ids_arr] = sizes_arr
    return dense_sizes


def parse_fasta_records(path):
    # Streams (header, color_set_id, seq_lines) from a '>unitig_id=N color_set_id=M ...'
    # Themisto2 FASTA export -- header/sequence lines are yielded unchanged so a
    # passing record can be written back out byte-for-byte.
    header = None
    color_set_id = None
    seq_lines = []
    with _open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, color_set_id, seq_lines
                header = line
                fields = dict(field.split("=") for field in line[1:].split())
                color_set_id = int(fields["color_set_id"])
                seq_lines = []
            else:
                seq_lines.append(line)
    if header is not None:
        yield header, color_set_id, seq_lines

##############################################################################
### Thresholding logic

def resolve_min_freq(raw: str) -> tuple[str, float]:
    # raw is a PRESET_MIN_FREQ name, or a custom fraction (e.g. '0.8') -- returns
    # (label, fraction); label is used as-is for presets, or turned filesystem-safe
    # ('0.8' -> 'freq0p8') for custom values, for output filenames.
    if raw in PRESET_MIN_FREQ:
        return raw, PRESET_MIN_FREQ[raw]
    try:
        freq = float(raw)
    except ValueError:
        raise ValueError(
            f"--min-freq must be 'core', 'relaxed', 'catchall', or a number "
            f"between 0.0 and 1.0 (got {raw!r})"
        )
    if not (0.0 <= freq <= 1.0):
        raise ValueError(f"--min-freq must be between 0.0 and 1.0, got {freq}")
    return f"freq{raw}".replace(".", "p"), freq


def compute_pass_masks(sizes: np.ndarray, n_colors: int, min_freq: float, min_genome_count: int) -> tuple[np.ndarray, np.ndarray]:
    # Two separate per-color_set_id boolean arrays (not one combined mask) so
    # filter_and_write_candidates() below can report WHY a unitig was dropped. A
    # unitig is a real candidate only where both are True.
    fractions = sizes / n_colors
    # >= min_freq, EXCEPT exactly at 0.0 (catchall's preset, or --min-freq 0): >= 0.0
    # is trivially true for everything, including a unitig completely absent from the
    # lineage (size == 0), which must never pass regardless of threshold.
    freq_pass = fractions > min_freq if min_freq == 0.0 else fractions >= min_freq
    count_pass = sizes >= min_genome_count
    return freq_pass, count_pass

##############################################################################
###  Filtering + output writing

def filter_and_write_candidates(unitigs_path, freq_pass, count_pass, out_path):
    # Single streaming pass: each record is written to out_path if its color_set_id
    # passes both masks, and tallied by drop reason otherwise. No intermediate ID
    # list, no second pass to re-extract sequences.
    total = kept = 0
    dropped_freq_only = dropped_count_only = dropped_both = 0
    with _open_out(out_path) as out_fh:
        for header, color_set_id, seq_lines in parse_fasta_records(unitigs_path):
            total += 1
            fp = color_set_id < len(freq_pass) and freq_pass[color_set_id]
            cp = color_set_id < len(count_pass) and count_pass[color_set_id]
            if fp and cp:
                kept += 1
                out_fh.write(header + "\n")
                out_fh.write("\n".join(seq_lines) + "\n")
            elif fp:
                dropped_count_only += 1
            elif cp:
                dropped_freq_only += 1
            else:
                dropped_both += 1

    return total, kept, dropped_freq_only, dropped_count_only, dropped_both


def write_lineage_summary_stats(lineage_id, min_freq_label, min_freq, min_genome_count,
                                 total, kept_count, dropped_freq_only, dropped_count_only,
                                 dropped_both, stats_out_path):
    dropped_total = dropped_freq_only + dropped_count_only + dropped_both
    stats_lines = [
        f"{'Lineage':<21} : {lineage_id}",
        f"{'Min freq':<21} : {min_freq_label} ({min_freq})",
        f"{'Min genome count':<21} : {min_genome_count}",
        f"{'Total unitigs scanned':<21} : {total:,}",
        f"{'Kept (E)':<21} : {kept_count:,}",
        f"{'Dropped (total)':<21} : {dropped_total:,}",
        f"{'  - fraction only':<21} : {dropped_freq_only:,}",
        f"{'  - count only':<21} : {dropped_count_only:,}",
        f"{'  - both':<21} : {dropped_both:,}",
    ]
    Path(stats_out_path).write_text("\n".join(stats_lines) + "\n")

##############################################################################
### Orchestrate

def main():
    args = parse_args()

    try:
        min_freq_label, min_freq = resolve_min_freq(args.min_freq)
    except ValueError as e:
        sys.exit(str(e))
    if args.min_genome_count < 0:
        sys.exit(f"--min-genome-count must be >= 0, got {args.min_genome_count}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.stats_output_dir:
        args.stats_output_dir.mkdir(parents=True, exist_ok=True)

    lineage_id = load_lineage_id(args.label_mapping)
    n_colors = load_n_colors(args.export_metadata)

    print(f"Loading color-set sizes from {args.color_sets} "
          f"(lineage={lineage_id}, n_colors={n_colors}) ...", file=sys.stderr)
    sizes = load_colorset_sizes(args.color_sets)
    print(f"  sizes loaded for {len(sizes):,} color_set_id(s)", file=sys.stderr)

    freq_pass, count_pass = compute_pass_masks(sizes, n_colors, min_freq, args.min_genome_count)

    out_path = args.output_dir / f"{lineage_id}_{min_freq_label}_candidate_unitigs.fasta"
    stats_out_path = (
        args.stats_output_dir / f"{lineage_id}_{min_freq_label}_stats.txt"
        if args.stats_output_dir else None
    )

    print(f"Filtering {args.unitigs} (lineage={lineage_id}, "
          f"min_freq={min_freq_label} ({min_freq}), "
          f"min_genome_count={args.min_genome_count}) ...", file=sys.stderr)
    total, kept_count, dropped_freq_only, dropped_count_only, dropped_both = filter_and_write_candidates(
        args.unitigs, freq_pass, count_pass, out_path
    )
    dropped_total = dropped_freq_only + dropped_count_only + dropped_both

    if stats_out_path:
        write_lineage_summary_stats(
            lineage_id, min_freq_label, min_freq, args.min_genome_count,
            total, kept_count, dropped_freq_only, dropped_count_only, dropped_both,
            stats_out_path,
        )

    print(f"  Total: {total:,}  Kept (E): {kept_count:,}  Dropped: {dropped_total:,} "
          f"(fraction only: {dropped_freq_only:,}, count only: {dropped_count_only:,}, "
          f"both: {dropped_both:,})", file=sys.stderr)
    if kept_count == 0:
        print(f"WARNING: no candidate unitigs survived filtering for lineage "
              f"'{lineage_id}' (min_freq={min_freq_label} ({min_freq}), "
              f"min_genome_count={args.min_genome_count}) -- {out_path} is empty.",
              file=sys.stderr)


if __name__ == "__main__":
    main()
