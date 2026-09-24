#!/usr/bin/env python3

"""
Filter a species-wide Themisto2 export down to a per-lineage candidate marker
unitig set -- keeping unitigs that are both lineage-CORE and lineage-SPECIFIC.

Supersedes core_catchall_filter.py (within-lineage core presets + genome-count
floor + FASTA output) and lineage_specificity_score.py (species-wide presence
scoring). Both are folded in here.

Runs on the SPECIES-wide export -- export.unitigs.fa / export.color_sets.txt plus
the species label_mapping.tsv (row order == colour ID order at build time) -- in
ONE streaming pass over color_sets.txt. That pass scores every lineage that can
affect a result: the requested --lineages plus every lineage large enough to
count as a sister in the outside max (>= --min-lineage-size). Singleton /
sub-threshold lineages are dropped from the presence matrix -- at species scale
(hundreds of one-genome clusters) that is the difference between a
(n_colour_sets x ~150) matrix and (n_colour_sets x every_label).

For each unitig / requested lineage:
  within_frac = (that lineage's genomes carrying the unitig) / (its genome count)
  outside     = the MAX, over every OTHER lineage with >= --min-lineage-size
                genomes, of THAT lineage's own presence fraction for the unitig
                (not a pooled "all other genomes" average -- a unitig concentrated
                in one sister lineage but absent everywhere else has a low pooled
                fraction but a high max fraction; pooling hides exactly the
                cross-reaction this filter exists to catch).

--- filter mode (default) -----------------------------------------------------
Keep a unitig for a lineage iff ALL of:
  * within_frac  >= --min-freq           ('core' 0.95 / 'relaxed' 0.5 /
                                           'catchall' >0 / a literal fraction)
  * lineage genome count >= --min-genome-count   (absolute floor, default 5)
  * --max-outside is unset  OR  outside <= --max-outside

--max-outside unset  ->  the outside constraint is not applied  ->  identical to
the old core-only filter. Set it (e.g. 0.05) to also require lineage-specificity.

Writes, per requested lineage, into --output-dir:
  <lineage>_candidate_unitigs.fasta   passing unitig sequences (empty + a warning
                                       if nothing passes -> caller skips the rebuild)
  <lineage>_specificity.tsv           unitig_id, within_pct, outside_pct, outside_lineage
  <lineage>_stats.txt                 kept / dropped counts + the thresholds used
                                       (in --stats-output-dir if given, else --output-dir)

--- score-only mode (--score-only) ------------------------------------------------
No filtering, no FASTA -- just <lineage>_specificity.tsv (the diagnostic view).

NB: the species graph and a lineage's own independently-built graph segment unitigs
differently, so this is not a literal re-derivation of the lineage-scoped filter's
kept/dropped set -- it's a self-consistent species-wide view, recompacted by the
downstream GGCAT/SBWT rebuild.
"""

import argparse
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

##############################################################################
# Presets / defaults

PRESET_MIN_FREQ = {"core": 0.95, "relaxed": 0.5, "catchall": 0.0}
DEFAULT_MIN_GENOME_COUNT = 5

##############################################################################
# CLI


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--unitigs", required=True, type=Path, help="Species-wide export.unitigs.fa")
    p.add_argument("--colour-sets", required=True, type=Path, help="Species-wide export.color_sets.txt")
    p.add_argument(
        "--label-mapping",
        required=True,
        type=Path,
        help="Species-wide label_mapping.tsv (Sample_ID, label -- ALL genomes; row order == colour ID)",
    )
    p.add_argument("--export-metadata", type=Path, default=None, help="Species-wide export.metadata.txt (num_colors)")
    p.add_argument(
        "--n-colours", type=int, default=None, help="Species genome count (alternative to --export-metadata)"
    )
    p.add_argument(
        "--lineages",
        nargs="*",
        default=[],
        help="Lineage label(s) to filter/score, e.g. --lineages 7PET. Mutually exclusive with --all-lineages.",
    )
    p.add_argument(
        "--all-lineages",
        action="store_true",
        help="Filter/score every lineage in --label-mapping that has >= --min-genome-count genomes "
        "(excluding 'unclassified'). Use when the caller requests no explicit --lineages.",
    )
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument(
        "--stats-output-dir", type=Path, default=None, help="Directory for <lineage>_stats.txt (default: --output-dir)"
    )
    p.add_argument(
        "-f",
        "--min-freq",
        default="core",
        help="Within-lineage presence-fraction cutoff: 'core' (>=0.95), 'relaxed' (>=0.5), "
        "'catchall' (>0), or a literal fraction e.g. 0.8. Default: core.",
    )
    p.add_argument(
        "--min-genome-count",
        type=int,
        default=DEFAULT_MIN_GENOME_COUNT,
        help=f"Absolute lineage-genome-count floor, on top of --min-freq (default {DEFAULT_MIN_GENOME_COUNT})",
    )
    p.add_argument(
        "--min-lineage-size",
        type=int,
        default=DEFAULT_MIN_GENOME_COUNT,
        help="Ignore sister lineages smaller than this when computing the outside max (default 5)",
    )
    p.add_argument(
        "--max-outside",
        type=float,
        default=None,
        help="Ceiling on the max-over-sister-lineage presence fraction. Unset -> not applied "
        "(core-only behaviour). Set e.g. 0.05 to also require lineage-specificity.",
    )
    p.add_argument(
        "--score-only",
        action="store_true",
        help="Diagnostic: write <lineage>_specificity.tsv only, no filtering, no FASTA.",
    )
    p.add_argument("--threads", type=int, default=4, help="Worker processes for the color_sets.txt streaming pass")
    return p.parse_args()


##############################################################################
# Loading


def _metadata_int(export_metadata_path, key):
    with open(export_metadata_path) as fh:
        for line in fh:
            k, _, value = line.strip().partition("=")
            if k == key:
                return int(value)
    return None


def load_n_colours(export_metadata_path) -> int:
    n = _metadata_int(export_metadata_path, "num_colors")
    if n is None:
        raise ValueError(f"{export_metadata_path} has no 'num_colors=' line")
    return n


def load_n_colour_sets(export_metadata_path):
    """num_colour_sets from export.metadata.txt (None if absent). Lets
    compute_presence_parallel pre-size `dense` instead of holding every chunk's
    counts until it can infer the row count."""
    return _metadata_int(export_metadata_path, "num_color_sets")


def load_colour_lineage(label_mapping_path, n_colours):
    """Returns (lineage_names, lineage_of_colour int array len n_colours, lineage_size
    int array). Row order of label_mapping == colour id (build-time convention)."""
    names, idx = [], {}
    lof = np.full(n_colours, -1, dtype=np.int64)
    with open(label_mapping_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        lc = header.index("label")
        for cid, line in enumerate(fh):
            lab = line.rstrip("\n").split("\t")[lc]
            if lab not in idx:
                idx[lab] = len(names)
                names.append(lab)
            if cid < n_colours:
                lof[cid] = idx[lab]
    sizes = np.bincount(lof[lof >= 0], minlength=len(names)).astype(np.int64)
    return names, lof, sizes


##############################################################################
# Streaming presence computation -- one pass over color_sets.txt, all lineages
# scored at once, chunked by byte offset across --threads worker processes.

_LOF = _NL = _PATH = None


def _init_worker(path, lof, nl):
    global _LOF, _NL, _PATH
    _PATH, _LOF, _NL = path, np.asarray(lof), nl


def _score_chunk(bounds):
    start, end = bounds
    csids, counts = [], []
    with open(_PATH, "rb") as f:
        f.seek(start)
        if start:
            f.readline()  # drop the partial line -- the previous chunk owns it
        while f.tell() < end:
            raw = f.readline()
            if not raw:
                break
            parts = raw.split()
            if len(parts) < 2:
                continue
            csids.append(int(parts[0].split(b"=")[1]))
            cids = np.array(parts[2:], dtype=np.int64)  # numpy C parser, not a Python int() comprehension
            # int32: a per-lineage count is bounded by that lineage's genome count
            # (<< 2**31). Halves the count matrix vs int64 -- the dominant cost of
            # this step at --all-lineages scale (n_colour_sets x hundreds).
            counts.append(np.bincount(_LOF[cids][_LOF[cids] >= 0], minlength=_NL).astype(np.int32))
    return np.array(csids, dtype=np.int64), (np.vstack(counts) if counts else np.zeros((0, _NL), dtype=np.int32))


def compute_presence_parallel(colour_sets_path, n_lineages, lof, threads, n_colour_sets=None):
    """Returns dense (n_colour_sets, n_lineages) int32 array: how many of each
    lineage's genomes are in each colour set.

    `n_lineages` is the RELEVANT lineage count (non-singletons + targets, see
    main), not every label -- that is the first memory lever. The second: when
    `n_colour_sets` is known (export.metadata.txt) `dense` is allocated up front
    and each worker chunk is scattered straight in, so peak is one `dense` plus a
    single chunk in flight, not `dense` plus the list of every chunk's counts."""
    size = colour_sets_path.stat().st_size
    n_chunks = max(threads * 3, 1)
    step = size // n_chunks
    bounds = [(i * step, size if i == n_chunks - 1 else (i + 1) * step) for i in range(n_chunks)]

    dense = np.zeros((n_colour_sets, n_lineages), dtype=np.int32) if n_colour_sets else None
    pending = []  # fallback path only (no metadata): (csid, counts) per chunk
    with ProcessPoolExecutor(
        max_workers=threads, initializer=_init_worker, initargs=(str(colour_sets_path), lof, n_lineages)
    ) as ex:
        for csid, cnt in ex.map(_score_chunk, bounds):
            if not len(csid):
                continue
            if dense is None:
                pending.append((csid, cnt))
                continue
            hi = int(csid.max()) + 1
            if hi > dense.shape[0]:  # metadata under-counted -- grow once
                dense = np.vstack([dense, np.zeros((hi - dense.shape[0], n_lineages), dtype=np.int32)])
            dense[csid] = cnt
    if dense is not None:
        return dense
    max_id = max((int(c.max()) for c, _ in pending), default=-1)
    dense = np.zeros((max_id + 1, n_lineages), dtype=np.int32)
    while pending:
        c, cnt = pending.pop()
        dense[c] = cnt
    return dense


##############################################################################
# Thresholding


def resolve_min_freq(raw: str):
    if raw in PRESET_MIN_FREQ:
        return raw, PRESET_MIN_FREQ[raw]
    try:
        freq = float(raw)
    except ValueError:
        raise ValueError(f"--min-freq must be 'core', 'relaxed', 'catchall', or a number 0.0-1.0 (got {raw!r})")
    if not (0.0 <= freq <= 1.0):
        raise ValueError(f"--min-freq must be between 0.0 and 1.0, got {freq}")
    return f"freq{raw}".replace(".", "p"), freq


# Columns of `dense` turned into a float fraction array at once when taking the
# outside max. Bounds the transient at (n_colour_sets x OUTSIDE_LINEAGE_BATCH)
# float32 regardless of how many sister lineages the species has.
OUTSIDE_LINEAGE_BATCH = 64


def lineage_view(dense, names, sizes, target_idx, min_lineage_size):
    """within_frac, within_count, max_outside_frac, max_outside_lineage_idx --
    all 1-D arrays indexed by colour_set_id, for one target lineage.

    The outside max is accumulated a column-batch at a time: the full
    (n_colour_sets x n_lineages) float copy this used to materialise -- once per
    target -- was the species-scale OOM."""
    n_cs, n_l = dense.shape
    within_count = dense[:, target_idx]
    with np.errstate(divide="ignore", invalid="ignore"):
        within_frac = np.nan_to_num(within_count / sizes[target_idx]) if sizes[target_idx] else np.zeros(n_cs)
        inv_size = np.where(sizes > 0, 1.0 / np.where(sizes > 0, sizes, 1), 0.0).astype(np.float32)

    sister = sizes >= min_lineage_size
    sister[target_idx] = False  # a lineage is never its own sister

    max_outside = np.zeros(n_cs, dtype=np.float32)
    max_outside_lineage = np.zeros(n_cs, dtype=np.int64)
    for s in range(0, n_l, OUTSIDE_LINEAGE_BATCH):
        cols = np.arange(s, min(s + OUTSIDE_LINEAGE_BATCH, n_l))
        cols = cols[sister[cols]]
        if not len(cols):
            continue
        fr = dense[:, cols].astype(np.float32)
        fr *= inv_size[cols]  # (n_cs, k), broadcast over colour sets
        b_max = fr.max(axis=1)
        b_arg = fr.argmax(axis=1)
        upd = b_max > max_outside
        max_outside[upd] = b_max[upd]
        max_outside_lineage[upd] = cols[b_arg[upd]]
    return within_frac, within_count, max_outside, max_outside_lineage


##############################################################################
# FASTA streaming


def parse_fasta_records(path):
    header = None
    colour_set_id = None
    seq_lines = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, colour_set_id, seq_lines
                header = line
                colour_set_id = int(dict(t.split("=") for t in line[1:].split())["color_set_id"])
                seq_lines = []
            else:
                seq_lines.append(line)
    if header is not None:
        yield header, colour_set_id, seq_lines


##############################################################################
# Output


def write_specificity_tsv(unitigs_path, within_frac, max_outside, max_outside_lineage, names, path):
    """Write specificity scores as percentages. Skip out-of-bounds unitigs."""
    n = len(within_frac)
    with open(path, "w") as tsv:
        tsv.write("unitig_id\twithin_pct\toutside_pct\toutside_lineage\n")
        for header, csid, _ in parse_fasta_records(unitigs_path):
            uid = header[1:].split()[0].split("=", 1)[-1]
            if csid < n:
                w = within_frac[csid] * 100
                o = max_outside[csid] * 100
                lin = names[int(max_outside_lineage[csid])]
                tsv.write(f"{uid}\t{w:.2f}\t{o:.2f}\t{lin}\n")


def filter_and_write(unitigs_path, keep, fasta_path, tsv_path, within_frac, max_outside, max_outside_lineage, names):
    """Filter unitigs, write FASTA + TSV with filtering decisions."""
    n = len(keep)
    total = kept = 0
    with open(fasta_path, "w") as fa, open(tsv_path, "w") as tsv:
        tsv.write("unitig_id\twithin_pct\toutside_pct\toutside_lineage\tkept\n")
        for header, csid, seq in parse_fasta_records(unitigs_path):
            total += 1
            uid = header[1:].split()[0].split("=", 1)[-1]
            ok = csid < n and keep[csid]
            w = (within_frac[csid] * 100) if csid < n else 0.0
            o = (max_outside[csid] * 100) if csid < n else 0.0
            lin = names[int(max_outside_lineage[csid])] if csid < n else "-"
            tsv.write(f"{uid}\t{w:.2f}\t{o:.2f}\t{lin}\t{int(bool(ok))}\n")
            if ok:
                kept += 1
                fa.write(header + "\n" + "\n".join(seq) + "\n")
    return total, kept


def write_stats(
    lineage_id, min_freq_label, min_freq, min_genome_count, max_outside_thresh, lineage_size, total, kept, path
):
    """Write filtering statistics."""
    lines = [
        f"{'Lineage':<24} : {lineage_id}",
        f"{'Lineage genome count':<24} : {lineage_size:,}",
        f"{'Min within-freq':<24} : {min_freq_label} ({min_freq})",
        f"{'Min genome count':<24} : {min_genome_count}",
        f"{'Max outside (max-over-lineage)':<24} : "
        f"{'not applied' if max_outside_thresh is None else max_outside_thresh}",
        f"{'Total unitigs scanned':<24} : {total:,}",
        f"{'Kept':<24} : {kept:,}",
    ]
    Path(path).write_text("\n".join(lines) + "\n")


##############################################################################
# Orchestrate


def main():
    args = parse_args()

    if args.export_metadata is None and args.n_colours is None:
        sys.exit("give --export-metadata or --n-colours")
    n_colours = args.n_colours if args.n_colours is not None else load_n_colours(args.export_metadata)
    n_colour_sets = load_n_colour_sets(args.export_metadata) if args.export_metadata is not None else None

    try:
        min_freq_label, min_freq = resolve_min_freq(args.min_freq)
    except ValueError as e:
        sys.exit(str(e))
    if args.min_genome_count < 0:
        sys.exit(f"--min-genome-count must be >= 0, got {args.min_genome_count}")
    if args.max_outside is not None and not (0.0 <= args.max_outside <= 1.0):
        sys.exit(f"--max-outside must be between 0.0 and 1.0, got {args.max_outside}")
    if bool(args.lineages) == bool(args.all_lineages):
        sys.exit("give exactly one of --lineages or --all-lineages")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    stats_dir = args.stats_output_dir or args.output_dir
    stats_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading lineage map from {args.label_mapping} ...", file=sys.stderr)
    names, lof, sizes = load_colour_lineage(args.label_mapping, n_colours)

    if args.all_lineages:
        targets = sorted(n for i, n in enumerate(names) if n != "unclassified" and sizes[i] >= args.min_genome_count)
        too_small = sorted(n for i, n in enumerate(names) if n != "unclassified" and sizes[i] < args.min_genome_count)
        if too_small:
            print(
                f"--all-lineages: skipping {len(too_small)} lineage(s) below --min-genome-count "
                f"({args.min_genome_count}): {', '.join(too_small)}",
                file=sys.stderr,
            )
        if not targets:
            print(
                f"--all-lineages: no lineage in label_mapping.tsv has >= {args.min_genome_count} genomes "
                "-- nothing to do.",
                file=sys.stderr,
            )
            return
    else:
        missing = [lid for lid in args.lineages if lid not in names]
        for lid in missing:
            print(f"WARNING: lineage '{lid}' not in label_mapping.tsv -- skipping", file=sys.stderr)
        targets = [lid for lid in args.lineages if lid in names]
        if not targets:
            sys.exit("none of the requested --lineages are in label_mapping.tsv -- nothing to do")

    # --- restrict the presence matrix to lineages that can actually matter ------
    # A column of `dense` is only ever read as (a) a requested target's own
    # within-count or (b) a sister lineage in the outside max -- and (b) already
    # ignores anything below --min-lineage-size. So keep only {size >=
    # min_lineage_size} plus the requested targets and remap colour->lineage
    # accordingly. For a species with hundreds of singleton lineages this is a
    # (n_colour_sets x ~150) matrix instead of (n_colour_sets x every_label).
    target_idx = [names.index(lid) for lid in targets]
    relevant = sizes >= args.min_lineage_size
    relevant[target_idx] = True
    keep_cols = np.flatnonzero(relevant)
    if len(keep_cols) < len(names):
        remap = np.full(len(names), -1, dtype=np.int64)
        remap[keep_cols] = np.arange(len(keep_cols))
        good = lof >= 0
        lof_r = np.full_like(lof, -1)
        lof_r[good] = remap[lof[good]]
        print(
            f"  presence matrix: {len(keep_cols)} relevant lineage(s) of {len(names)} "
            f"(dropped {len(names) - len(keep_cols)} below --min-lineage-size {args.min_lineage_size})",
            file=sys.stderr,
        )
        names = [names[i] for i in keep_cols]
        sizes = sizes[keep_cols]
        lof = lof_r

    print(
        f"Streaming {args.colour_sets} across {args.threads} worker(s), scoring {len(targets)} lineage(s) ...",
        file=sys.stderr,
    )
    dense = compute_presence_parallel(args.colour_sets, len(names), lof, args.threads, n_colour_sets=n_colour_sets)

    for lid in targets:
        t = names.index(lid)
        within_frac, within_count, max_outside, max_outside_lineage = lineage_view(
            dense, names, sizes, t, args.min_lineage_size
        )

        if args.score_only:
            tsv = args.output_dir / f"{lid}_specificity.tsv"
            write_specificity_tsv(args.unitigs, within_frac, max_outside, max_outside_lineage, names, tsv)
            print(f"  wrote {tsv}", file=sys.stderr)
            continue

        freq_pass = (within_frac > 0.0) if min_freq == 0.0 else (within_frac >= min_freq)
        count_pass = within_count >= args.min_genome_count
        outside_pass = (
            np.ones_like(within_frac, dtype=np.bool_) if args.max_outside is None else (max_outside <= args.max_outside)
        )
        keep = freq_pass & count_pass & outside_pass

        fasta = args.output_dir / f"{lid}_candidate_unitigs.fasta"
        tsv = args.output_dir / f"{lid}_specificity.tsv"
        total, kept = filter_and_write(
            args.unitigs, keep, fasta, tsv, within_frac, max_outside, max_outside_lineage, names
        )
        write_stats(
            lid,
            min_freq_label,
            min_freq,
            args.min_genome_count,
            args.max_outside,
            int(sizes[t]),
            total,
            kept,
            stats_dir / f"{lid}_stats.txt",
        )
        print(f"  {lid}: scanned {total:,}  kept {kept:,}", file=sys.stderr)
        if kept == 0:
            print(f"WARNING: no candidate unitigs survived for lineage '{lid}' -- {fasta} is empty.", file=sys.stderr)


if __name__ == "__main__":
    main()