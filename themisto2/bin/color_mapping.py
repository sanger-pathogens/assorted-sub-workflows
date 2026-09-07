#!/usr/bin/env python3
"""
Build Themisto2's ``--file-colors`` input for a species-wide colour index.

Reads one metadata table (one row per genome) and writes, into --output_dir,
three files prefixed with the species name (--species-name):

  <species>_file_colors_input.txt  assembly paths, one per line, grouped by label
  <species>_label_mapping.tsv      Sample_ID -> label, in index (colour-ID) order
  <species>_stats.json             summary counts, incl. assemblies per group
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from pathlib import Path

import pandas as pd

UNCLASSIFIED = "unclassified"

# Sanitisation applied when the assembly FASTA files were written to disk.
_UNSAFE = re.compile(r"[^A-Za-z0-9._-]")


def fs_safe(name: str) -> str:
    """Filesystem-safe form used to match a metadata Sample_ID to a FASTA name."""
    return _UNSAFE.sub("_", name)


def parse_args():
    p = argparse.ArgumentParser(description="Build Themisto2 file-colors input from a metadata table.")
    p.add_argument("--metadata", required=True, help="Path to the metadata CSV.")
    p.add_argument(
        "--species-name",
        required=True,
        help="Species identifier (from the manifest); used as the output filename prefix.",
    )
    p.add_argument(
        "--sample-col",
        default="Sample_ID",
        help="Metadata column holding sample identifiers (default: Sample_ID).",
    )
    p.add_argument(
        "--group-label",
        required=True,
        help="Metadata column to group and sort genomes by, written to label_mapping.tsv.",
    )
    asm = p.add_mutually_exclusive_group(required=True)
    asm.add_argument("--assembly-dir", help="Directory containing assembly FASTA files.")
    asm.add_argument("--assembly-paths", help="Text file listing one assembly path per line.")
    p.add_argument(
        "--assembly-suffix",
        default=".contigs.fasta",
        help="Suffix appended to Sample_ID to form the assembly filename (default: .contigs.fasta).",
    )
    p.add_argument("--output_dir", required=True, help="Directory to write the output files into.")
    return p.parse_args()


def load_assemblies(assembly_input: str) -> tuple[dict, list[str]]:
    """(basename -> full path, dead paths). From a directory or a path-list file.

    Path-list entries that don't point at an existing file are pulled out into
    the second list -- a wrong path in the .txt must not reach file_colors_input.
    """
    path = Path(assembly_input)
    if path.is_dir():
        # Resolve to a canonical absolute path so file_colors_input.txt is readable
        # from a different task's work dir (staged symlinks won't resolve there).
        root = path.resolve()
        return {name: str(root / name) for name in os.listdir(root)}, []
    assemblies: dict = {}
    dead: list[str] = []
    with open(path) as fh:
        for line in fh:
            entry = line.strip()
            if not entry:
                continue
            if Path(entry).is_file():
                assemblies[Path(entry).name] = entry
            else:
                dead.append(entry)
    return assemblies, dead


def index_by_safe_name(assemblies: dict) -> dict:
    """basename->path keyed instead by fs_safe(basename); first wins on a clash."""
    safe: dict = {}
    for basename, full_path in assemblies.items():
        key = fs_safe(basename)
        if key in safe:
            print(
                f"warning: '{basename}' and '{Path(safe[key]).name}' normalise to the "
                f"same name -- keeping the first",
                file=sys.stderr,
            )
            continue
        safe[key] = full_path
    return safe


def strip_suffix(basename: str, suffix: str) -> str:
    return basename[: -len(suffix)] if suffix and basename.endswith(suffix) else Path(basename).stem


def main():
    args = parse_args()
    prefix = args.species_name

    metadata = pd.read_csv(args.metadata, low_memory=False)

    missing = [c for c in (args.sample_col, args.group_label) if c not in metadata.columns]
    if missing:
        sys.exit(f"Error: column(s) not found in metadata: {', '.join(missing)}")

    metadata[args.sample_col] = metadata[args.sample_col].astype(str)
    duplicates = int(metadata[args.sample_col].duplicated().sum())
    if duplicates:
        metadata = metadata.drop_duplicates(subset=args.sample_col)
        print(f"warning: ignored {duplicates} row(s) with a duplicate {args.sample_col}", file=sys.stderr)

    # Genomes with no group label aren't dropped -- relabel and keep them.
    blank_label = metadata[args.group_label].isna()
    metadata[args.group_label] = metadata[args.group_label].fillna(UNCLASSIFIED).astype(str)

    # Match each metadata row to an assembly by normalised filename.
    assemblies, dead_paths = load_assemblies(args.assembly_dir or args.assembly_paths)
    if dead_paths:
        print(
            f"warning: {len(dead_paths)} assembly path(s) don't point at a file, e.g. "
            + ", ".join(dead_paths[:3]),
            file=sys.stderr,
        )
    by_safe_name = index_by_safe_name(assemblies)
    metadata["_want"] = (metadata[args.sample_col] + args.assembly_suffix).map(fs_safe)
    matched = metadata["_want"].isin(by_safe_name)

    kept = metadata[matched].copy()
    kept["file_path"] = kept["_want"].map(by_safe_name)
    n_dropped = int((~matched).sum())
    n_relabelled = int((blank_label & matched).sum())

    # Assembly files no metadata row claimed -> keep them too, as "unclassified".
    claimed = set(kept["_want"])
    orphan_paths = sorted(p for safe, p in by_safe_name.items() if safe not in claimed)
    orphan_rows = pd.DataFrame(
        {
            args.sample_col: [strip_suffix(Path(p).name, args.assembly_suffix) for p in orphan_paths],
            args.group_label: UNCLASSIFIED,
            "file_path": orphan_paths,
        }
    )

    written = pd.concat(
        [kept[[args.sample_col, args.group_label, "file_path"]], orphan_rows], ignore_index=True
    )
    # sort by label, then sample ID within each label -- fixes colour-ID order
    written = written.sort_values(
        [args.group_label, args.sample_col], key=lambda c: c.astype(str)
    )

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    written["file_path"].to_csv(out / f"{prefix}_file_colors_input.txt", index=False, header=False)
    (
        written[[args.sample_col, args.group_label]]
        .rename(columns={args.sample_col: "Sample_ID", args.group_label: "label"})
        .to_csv(out / f"{prefix}_label_mapping.tsv", index=False, sep="\t")
    )

    per_group = {
        str(label): int(n)
        for label, n in written[args.group_label].value_counts().sort_index().items()
    }
    stats = {
        "species": prefix,
        "group_label_column": args.group_label,
        "assemblies_total": len(written) + n_dropped,
        "assemblies_written": len(written),
        "assemblies_dropped_fasta_not_found": n_dropped,
        "assembly_paths_missing_file": len(dead_paths),
        "relabelled_unclassified": n_relabelled,
        "assemblies_without_metadata_row": len(orphan_paths),
        "assemblies_per_group": per_group,
        "note": (
            "'unclassified' = metadata rows with a blank group label plus assembly "
            "files with no metadata row; both are kept in the index. "
            "assemblies_dropped_fasta_not_found are the only genomes left out."
        ),
    }
    with open(out / f"{prefix}_stats.json", "w") as fh:
        json.dump(stats, fh, indent=2)

    summary = f"{prefix}: {stats['assemblies_written']}/{stats['assemblies_total']} assemblies written"
    if n_dropped:
        summary += f"; {n_dropped} dropped (FASTA not found)"
    if orphan_paths:
        summary += f"; {len(orphan_paths)} had no metadata row -> {UNCLASSIFIED}"
    if n_relabelled:
        summary += f"; {n_relabelled} blank label -> {UNCLASSIFIED}"
    print(summary, file=sys.stderr)


if __name__ == "__main__":
    main()
