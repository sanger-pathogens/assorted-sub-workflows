#!/usr/bin/env python3

"""
Scores unitigs by within-lineage vs outside-lineage presence, using the
SPECIES-wide Themisto2 export (all lineages together, not the
lineage-scoped exports CANDIDATE_FILTER uses) -- a diagnostic view of
specificity, meant to be looked at before committing to expensive
index-native set-diff operations, not a filtering step itself.

within_pct/outside_pct for a color_set_id are computed via a lineage
membership mask over the species-wide colour space (species-wide
label_mapping.tsv's row order == colour ID order at build time), batched
across every --lineages value in one streaming pass over
export.color_sets.txt. specificity_score = within_pct - outside_pct
(Youden's J).

Caveat: the species-wide graph and each lineage's own independently-built
graph (CANDIDATE_FILTER's input) can segment unitigs differently, so this
is NOT a literal re-derivation of CANDIDATE_FILTER's own kept/dropped
decision for that lineage (different build, different unitig IDs) -- it's
a separate, self-consistent species-wide view.

Output: '{lineage_id}_specificity.tsv' per lineage, plus
'{lineage_id}_specificity.png' if --plot (skipped, not an error, if
matplotlib isn't installed).
"""

import argparse
import gzip
import sys
from pathlib import Path
import numpy as np

##############################################################################
### CLI

def parse_args():
    parser = argparse.ArgumentParser(
        description="Score unitigs by within-lineage vs outside-lineage "
                     "presence using the SPECIES-wide Themisto2 export, for "
                     "one or more lineages at once -- a diagnostic view, "
                     "independent of (and not gating) CANDIDATE_FILTER's "
                     "own lineage-scoped filtering.",
    )
    parser.add_argument("--unitigs", required=True, type=Path,
                         help="Species-wide export.unitigs.fa(.gz).")
    parser.add_argument("--color-sets", required=True, type=Path,
                         help="Species-wide export.color_sets.txt(.gz).")
    parser.add_argument("--label-mapping", required=True, type=Path,
                         help="Species-wide label_mapping.tsv (Sample_ID, "
                              "label -- covering ALL genomes, not one "
                              "lineage's subset).")
    parser.add_argument("--n-colors", required=True, type=int,
                         help="Species-wide genome count (species "
                              "export.metadata.txt's num_colors).")
    parser.add_argument("--lineages", required=True, nargs="+",
                         help="Lineage label(s) to score, e.g. --lineages "
                              "GPSC1 GPSC2 -- matches --target_groups.")
    parser.add_argument("--output-dir", required=True, type=Path,
                         help="Directory for '{lineage_id}_specificity.tsv' "
                              "(+ '{lineage_id}_specificity.png' if --plot).")
    parser.add_argument("--plot", action="store_true",
                         help="Also generate a specificity plot per lineage "
                              "(requires matplotlib; skipped with a warning, "
                              "not an error, if it's not installed).")
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
### I/O loading

def load_lineage_map(label_mapping_path) -> dict[str, list[int]]:
    # label_mapping.tsv's row order == colour ID order at build time (0-indexed,
    # header on line 1) -- same convention color_mapping.py itself relies on.
    lineage_map: dict[str, list[int]] = {}
    with _open(label_mapping_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        label_col = header.index("label")
        for color_id, line in enumerate(fh):
            fields = line.rstrip("\n").split("\t")
            lineage_map.setdefault(fields[label_col], []).append(color_id)
    return lineage_map


def load_unitig_colorset_ids(unitigs_path) -> tuple[list[str], np.ndarray]:
    # Sequences aren't needed for scoring -- only each unitig's own ID (for
    # labeling rows) and its color_set_id.
    unitig_ids = []
    colorset_ids = []
    with _open(unitigs_path) as f:
        for line in f:
            if line.startswith(">"):
                fields = dict(field.split("=") for field in line[1:].split())
                unitig_ids.append(fields["unitig_id"])
                colorset_ids.append(int(fields["color_set_id"]))
    return unitig_ids, np.array(colorset_ids, dtype=np.int64)

##############################################################################
### Fraction computation

def build_lineage_masks(lineage_map, lineages, n_colors) -> dict[str, np.ndarray]:
    masks = {}
    for lineage_id in lineages:
        if lineage_id not in lineage_map:
            print(f"WARNING: lineage '{lineage_id}' not found in label_mapping.tsv -- skipping",
                  file=sys.stderr)
            continue
        mask = np.zeros(n_colors, dtype=np.bool_)
        mask[lineage_map[lineage_id]] = True
        masks[lineage_id] = mask
    return masks


def compute_fractions_streaming(color_sets_path, lineage_masks: dict[str, np.ndarray]) -> tuple[dict, dict]:
    # Single streaming pass over the (potentially huge) species-wide
    # color_sets.txt, computing within/outside presence for every requested
    # lineage at once -- batched so the file is only ever read once, not once
    # per lineage.
    lineage_sizes = {lid: int(np.sum(m)) for lid, m in lineage_masks.items()}
    n_total = next(iter(lineage_masks.values())).shape[0]
    outside_masks = {lid: ~m for lid, m in lineage_masks.items()}
    outside_sizes = {lid: n_total - lineage_sizes[lid] for lid in lineage_masks}

    seen_ids = []
    within_by_lineage = {lid: [] for lid in lineage_masks}
    outside_by_lineage = {lid: [] for lid in lineage_masks}

    with _open(color_sets_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            color_set_id = int(parts[0].split("=")[1])
            color_ids = np.array([int(tok) for tok in parts[2:]], dtype=np.int64)
            seen_ids.append(color_set_id)
            for lid, mask in lineage_masks.items():
                within_by_lineage[lid].append(np.sum(mask[color_ids]) / lineage_sizes[lid])
                outside_by_lineage[lid].append(np.sum(outside_masks[lid][color_ids]) / outside_sizes[lid])

    ids_arr = np.array(seen_ids, dtype=np.int64)
    max_id = int(ids_arr.max()) if len(ids_arr) else -1

    def to_dense(by_lineage):
        dense_by_lineage = {}
        for lid, values in by_lineage.items():
            dense = np.zeros(max_id + 1, dtype=np.float64)
            dense[ids_arr] = np.array(values, dtype=np.float64)
            dense_by_lineage[lid] = dense
        return dense_by_lineage

    return to_dense(within_by_lineage), to_dense(outside_by_lineage)

##############################################################################
### Output writing

def write_specificity_tsv(unitig_ids, colorset_ids, within_fractions, outside_fractions, out_path):
    with _open_out(out_path) as out_fh:
        out_fh.write("unitig_id\twithin_pct\toutside_pct\tspecificity_score\n")
        for uid, csid in zip(unitig_ids, colorset_ids):
            w = within_fractions[csid] if csid < len(within_fractions) else 0.0
            o = outside_fractions[csid] if csid < len(outside_fractions) else 0.0
            out_fh.write(f"{uid}\t{w:.6f}\t{o:.6f}\t{w - o:.6f}\n")


def plot_specificity(within_fractions, outside_fractions, colorset_ids, lineage_id, out_path):
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("WARNING: matplotlib not installed -- skipping --plot.", file=sys.stderr)
        return

    within = within_fractions[colorset_ids]
    outside = outside_fractions[colorset_ids]
    score = within - outside

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(f"Lineage {lineage_id} specificity ({len(colorset_ids):,} unitigs)",
                 fontsize=13, fontweight="bold")

    # Hexbin, not a raw scatter -- unitig counts can run into the hundreds of
    # thousands, where a scatter of that many points is both slow and unreadable.
    ax = axes[0]
    hb = ax.hexbin(outside, within, gridsize=60, cmap="viridis", mincnt=1)
    ax.plot([0, 1], [0, 1], color="orange", linestyle="--", linewidth=1, label="within == outside")
    ax.set_xlabel("Outside-lineage presence (outside_pct)")
    ax.set_ylabel("Within-lineage presence (within_pct)")
    ax.set_title("Presence: within vs outside")
    ax.legend()
    plt.colorbar(hb, ax=ax, label="unitig count")

    ax = axes[1]
    ax.hist(score, bins=50, color="steelblue", edgecolor="black", alpha=0.7)
    ax.axvline(0, color="orange", linestyle="--", linewidth=1)
    ax.set_xlabel("Specificity score (within_pct - outside_pct)")
    ax.set_ylabel("Unitig count")
    ax.set_title("Specificity score distribution")

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    print(f"Saved plot to {out_path}", file=sys.stderr)

##############################################################################
### Orchestrate

def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading lineage map from {args.label_mapping} ...", file=sys.stderr)
    lineage_map = load_lineage_map(args.label_mapping)

    lineage_masks = build_lineage_masks(lineage_map, args.lineages, args.n_colors)
    if not lineage_masks:
        sys.exit("No requested lineage found in label_mapping.tsv -- nothing to score.")

    print(f"Streaming color sets from {args.color_sets} for {len(lineage_masks)} lineage(s) ...",
          file=sys.stderr)
    within_by_lineage, outside_by_lineage = compute_fractions_streaming(args.color_sets, lineage_masks)

    print(f"Loading unitig color_set_ids from {args.unitigs} ...", file=sys.stderr)
    unitig_ids, colorset_ids = load_unitig_colorset_ids(args.unitigs)
    print(f"  {len(unitig_ids):,} unitig(s) loaded", file=sys.stderr)

    for lineage_id in lineage_masks:
        out_path = args.output_dir / f"{lineage_id}_specificity.tsv"
        write_specificity_tsv(
            unitig_ids, colorset_ids, within_by_lineage[lineage_id], outside_by_lineage[lineage_id], out_path
        )
        print(f"  Wrote {out_path}", file=sys.stderr)

        if args.plot:
            plot_path = args.output_dir / f"{lineage_id}_specificity.png"
            plot_specificity(
                within_by_lineage[lineage_id], outside_by_lineage[lineage_id], colorset_ids, lineage_id, plot_path
            )


if __name__ == "__main__":
    main()
