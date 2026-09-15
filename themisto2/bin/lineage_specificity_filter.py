#!/usr/bin/env python3

"""
Filter a species-wide Themisto2 export down to a per-lineage candidate marker
unitig set -- keeping unitigs that are both lineage-CORE and lineage-SPECIFIC.

Supersedes core_catchall_filter.py (within-lineage core presets + genome-count
floor + FASTA output) and lineage_specificity_score.py (species-wide presence
scoring). Both are folded in here.

Runs on the SPECIES-wide export -- export.unitigs.fa / export.color_sets.txt plus
the species label_mapping.tsv (row order == colour ID order at build time) -- in
ONE streaming pass over color_sets.txt that scores every lineage present, however
many --lineages are actually requested.

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
    p.add_argument("--color-sets", required=True, type=Path, help="Species-wide export.color_sets.txt")
    p.add_argument(
        "--label-mapping",
        required=True,
        type=Path,
        help="Species-wide label_mapping.tsv (Sample_ID, label -- ALL genomes; row order == colour ID)",
    )
    p.add_argument("--export-metadata", type=Path, default=None, help="Species-wide export.metadata.txt (num_colors)")
    p.add_argument("--n-colors", type=int, default=None, help="Species genome count (alternative to --export-metadata)")
    p.add_argument(
        "--lineages", required=True, nargs="+", help="Lineage label(s) to filter/score, e.g. --lineages 7PET"
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


def load_n_colors(export_metadata_path) -> int:
    with open(export_metadata_path) as fh:
        for line in fh:
            key, _, value = line.strip().partition("=")
            if key == "num_colors":
                return int(value)
    raise ValueError(f"{export_metadata_path} has no 'num_colors=' line")


def load_colour_lineage(label_mapping_path, n_colors):
    """Returns (lineage_names, lineage_of_colour int array len n_colors, lineage_size
    int array). Row order of label_mapping == colour id (build-time convention)."""
    names, idx = [], {}
    lof = np.full(n_colors, -1, dtype=np.int64)
    with open(label_mapping_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        lc = header.index("label")
        for cid, line in enumerate(fh):
            lab = line.rstrip("\n").split("\t")[lc]
            if lab not in idx:
                idx[lab] = len(names)
                names.append(lab)
            if cid < n_colors:
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
            counts.append(np.bincount(_LOF[cids][_LOF[cids] >= 0], minlength=_NL).astype(np.int64))
    return np.array(csids, dtype=np.int64), (np.vstack(counts) if counts else np.zeros((0, _NL), dtype=np.int64))


def compute_presence_parallel(color_sets_path, n_lineages, lof, threads):
    """Returns dense (n_color_sets, n_lineages) int array: how many of each
    lineage's genomes are in each colour set."""
    size = color_sets_path.stat().st_size
    n_chunks = max(threads * 3, 1)
    step = size // n_chunks
    bounds = [(i * step, size if i == n_chunks - 1 else (i + 1) * step) for i in range(n_chunks)]

    all_csid, all_counts = [], []
    with ProcessPoolExecutor(
        max_workers=threads, initializer=_init_worker, initargs=(str(color_sets_path), lof, n_lineages)
    ) as ex:
        for csid, cnt in ex.map(_score_chunk, bounds):
            all_csid.append(csid)
            all_counts.append(cnt)
    csid = np.concatenate(all_csid) if all_csid else np.zeros(0, dtype=np.int64)
    max_id = int(csid.max()) if len(csid) else -1
    dense = np.zeros((max_id + 1, n_lineages), dtype=np.int64)
    dense[csid] = np.vstack(all_counts) if all_counts else np.zeros((0, n_lineages), dtype=np.int64)
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


def lineage_view(dense, names, sizes, target_idx, min_lineage_size):
    """within_frac, within_count, max_outside_frac, max_outside_lineage_idx --
    all dense arrays indexed by color_set_id, for one target lineage."""
    within_count = dense[:, target_idx]
    with np.errstate(divide="ignore", invalid="ignore"):
        within_frac = np.nan_to_num(within_count / sizes[target_idx]) if sizes[target_idx] else np.zeros(len(dense))
        fr = np.nan_to_num(dense / np.where(sizes > 0, sizes, 1))
    fr = fr.copy()
    fr[:, target_idx] = 0.0
    fr[:, sizes < min_lineage_size] = 0.0
    max_outside = fr.max(axis=1)
    max_outside_lineage = fr.argmax(axis=1)
    return within_frac, within_count, max_outside, max_outside_lineage


##############################################################################
# FASTA streaming


def parse_fasta_records(path):
    header = None
    color_set_id = None
    seq_lines = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, color_set_id, seq_lines
                header = line
                color_set_id = int(dict(t.split("=") for t in line[1:].split())["color_set_id"])
                seq_lines = []
            else:
                seq_lines.append(line)
    if header is not None:
        yield header, color_set_id, seq_lines


##############################################################################
# Output


def write_specificity_tsv(unitigs_path, within_frac, max_outside, max_outside_lineage, names, path):
    n = len(within_frac)
    with open(path, "w") as tsv:
        tsv.write("unitig_id\twithin_pct\toutside_pct\toutside_lineage\n")
        for header, csid, _ in parse_fasta_records(unitigs_path):
            uid = header[1:].split()[0].split("=", 1)[-1]
            if csid < n:
                tsv.write(
                    f"{uid}\t{within_frac[csid]:.6f}\t{max_outside[csid]:.6f}\t{names[max_outside_lineage[csid]]}\n"
                )
            else:
                tsv.write(f"{uid}\t0.000000\t0.000000\t-\n")


def filter_and_write(unitigs_path, keep, fasta_path, tsv_path, within_frac, max_outside, max_outside_lineage, names):
    n = len(keep)
    total = kept = 0
    with open(fasta_path, "w") as fa, open(tsv_path, "w") as tsv:
        tsv.write("unitig_id\twithin_pct\toutside_pct\toutside_lineage\tkept\n")
        for header, csid, seq in parse_fasta_records(unitigs_path):
            total += 1
            uid = header[1:].split()[0].split("=", 1)[-1]
            ok = csid < n and keep[csid]
            w = within_frac[csid] if csid < n else 0.0
            o = max_outside[csid] if csid < n else 0.0
            lin = names[max_outside_lineage[csid]] if csid < n else "-"
            tsv.write(f"{uid}\t{w:.6f}\t{o:.6f}\t{lin}\t{int(bool(ok))}\n")
            if ok:
                kept += 1
                fa.write(header + "\n" + "\n".join(seq) + "\n")
    return total, kept


def write_stats(
    lineage_id, min_freq_label, min_freq, min_genome_count, max_outside_thresh, lineage_size, total, kept, path
):
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

    if args.export_metadata is None and args.n_colors is None:
        sys.exit("give --export-metadata or --n-colors")
    n_colors = args.n_colors if args.n_colors is not None else load_n_colors(args.export_metadata)

    try:
        min_freq_label, min_freq = resolve_min_freq(args.min_freq)
    except ValueError as e:
        sys.exit(str(e))
    if args.min_genome_count < 0:
        sys.exit(f"--min-genome-count must be >= 0, got {args.min_genome_count}")
    if args.max_outside is not None and not (0.0 <= args.max_outside <= 1.0):
        sys.exit(f"--max-outside must be between 0.0 and 1.0, got {args.max_outside}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    stats_dir = args.stats_output_dir or args.output_dir
    stats_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading lineage map from {args.label_mapping} ...", file=sys.stderr)
    names, lof, sizes = load_colour_lineage(args.label_mapping, n_colors)
    missing = [lid for lid in args.lineages if lid not in names]
    if missing:
        for lid in missing:
            print(f"WARNING: lineage '{lid}' not in label_mapping.tsv -- skipping", file=sys.stderr)
    targets = [lid for lid in args.lineages if lid in names]
    if not targets:
        sys.exit("none of the requested --lineages are in label_mapping.tsv -- nothing to do")

    print(
        f"Streaming {args.color_sets} across {args.threads} worker(s), scoring {len(targets)} lineage(s) ...",
        file=sys.stderr,
    )
    dense = compute_presence_parallel(args.color_sets, len(names), lof, args.threads)

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
