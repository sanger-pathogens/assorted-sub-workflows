#!/usr/bin/env python3
"""
Generic color mapping script for Themisto2 --file-colors input.

Reads metadata CSV, sorts by group label column, matches samples to
assemblies, and prepares index input files organized in subdirectories:
  - index_species/ : whole-species index files
  - index_target_group/<group>/ : per-group indexes (if --target-groups provided)

The metadata layout is not fixed: any CSV/TSV works as long as it has one
row per genome, a sample-identifier column (--sample-col) and a grouping
column (--group-label). The grouping column can be any level -- lineage,
sub-lineage, strain, serotype -- the script does not interpret it.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Build Themisto2 file-colors input from a metadata table.")
    parser.add_argument("--metadata", required=True, help="Path to metadata CSV file.")
    parser.add_argument(
        "--sample-col",
        default="Sample_ID",
        help="Column name in metadata containing sample identifiers (default: Sample_ID).",
    )
    parser.add_argument(
        "--group-label",
        required=True,
        help="Column name to sort samples by (lexicographically), write to "
             "label_mapping.tsv, and match against --target-groups.",
    )
    parser.add_argument(
        "--target-groups",
        default="",
        help=(
            "Comma-separated label(s) to build separate lineage-scoped indexes for, "
            "e.g. GPSC1,GPSC2. Leave empty for species-wide only."
        ),
    )
    asm = parser.add_mutually_exclusive_group(required=True)
    asm.add_argument("--assembly-dir", help="Directory containing assembly FASTA files.")
    asm.add_argument("--assembly-paths", help="Text file listing one assembly path per line.")
    parser.add_argument(
        "--assembly-suffix",
        default=".contigs.fasta",
        help=(
            "Suffix appended to Sample_ID to form the assembly filename "
            "(default: .contigs.fasta). Only used with --assembly-dir."
        ),
    )
    parser.add_argument("--output_dir", required=True, help="Output directory for index subdirectories.")
    return parser.parse_args()


def match_assemblies(metadata, assembly_input, assembly_suffix, sample_col, group_label):
    """Match metadata samples to assembly files on disk."""
    assembly_input_path = Path(assembly_input)

    if assembly_input_path.is_dir():
        # Resolve to canonical absolute path so file_colors_input.txt can be read
        # by GGCAT in a different task directory (symlinks won't resolve there)
        assembly_dir = assembly_input_path.resolve()
        assembly_files = set(os.listdir(assembly_dir))
        metadata["filename"] = metadata[sample_col].astype(str) + assembly_suffix
        no_assembly = ~metadata["filename"].isin(assembly_files)
        no_assembly_by_group = metadata.loc[no_assembly, group_label].astype(str).value_counts().to_dict()
        on_disk_no_metadata = assembly_files - set(metadata["filename"])
        metadata = metadata[metadata["filename"].isin(assembly_files)].reset_index(drop=True)
        metadata["file_path"] = metadata["filename"].apply(lambda f: str(assembly_dir / f))
    else:
        with open(assembly_input_path) as fh:
            path_list = [line.strip() for line in fh if line.strip()]
        # build lookup: basename (without suffix) -> full path
        path_lookup = {Path(p).name: p for p in path_list}
        metadata["filename"] = metadata[sample_col].astype(str) + assembly_suffix
        no_assembly = ~metadata["filename"].isin(path_lookup)
        no_assembly_by_group = metadata.loc[no_assembly, group_label].astype(str).value_counts().to_dict()
        on_disk_no_metadata = set(path_lookup.keys()) - set(metadata["filename"])
        metadata = metadata[metadata["filename"].isin(path_lookup)].reset_index(drop=True)
        metadata["file_path"] = metadata["filename"].apply(lambda f: path_lookup[f])

    return metadata, no_assembly_by_group, on_disk_no_metadata


def write_index_files(
    metadata: pd.DataFrame,
    output_dir: Path,
    index_type: str,
    group_name: str | None,
    group_label: str,
    sample_col: str,
    nan_label: int | None = None,
    no_assembly: int | None = None,
    on_disk_no_metadata: int | None = None,
) -> None:
    """Write file_colors_input.txt, label_mapping.tsv and stats.json to the appropriate subdirectory.

    Args:
        metadata: Rows for this index only (already filtered/sorted), with a
            ``file_path`` column plus ``sample_col`` and ``group_label``.
        output_dir: Root output dir; files go under ``index_species/`` or
            ``index_target_group/<group_name>/``.
        index_type: ``"species"`` for the whole-species index, anything else
            for a per-group index.
        group_name: Target-group label for a per-group index; ``None`` for the
            species index.
        group_label: Metadata column used for grouping/labelling.
        sample_col: Metadata column holding sample identifiers.
        nan_label: Count of samples dropped for a missing group label. Only
            meaningful species-wide; pass ``None`` for a per-group index so it
            is omitted from stats.json (not attributable to one group).
        no_assembly: Count of samples dropped for having no assembly on disk.
        on_disk_no_metadata: Count of assemblies on disk with no metadata row.
            Species-wide only; ``None`` for a per-group index.
    """
    if index_type == "species":
        index_dir = output_dir / "index_species"
        tag = "species"
    else:
        if not group_name or group_name in (".", "..") or "/" in group_name or "\\" in group_name:
            sys.exit(f"Error: invalid --target-groups entry {group_name!r} -- must not contain path separators.")
        index_dir = output_dir / "index_target_group" / group_name
        tag = group_name

    index_dir.mkdir(parents=True, exist_ok=True)

    metadata["file_path"].to_csv(index_dir / f"{tag}_file_colors_input.txt", index=False, header=False)

    label_mapping = (
        metadata[[sample_col, group_label]]
        .drop_duplicates()
        .rename(columns={sample_col: "Sample_ID", group_label: "label"})
    )
    label_mapping.to_csv(index_dir / f"{tag}_label_mapping.tsv", index=False, sep="\t")

    stats = {
        "index_type": index_type,
        "target_group": group_name or "species_wide",
        "metadata_column": group_label,
    }
    # nan_label/on_disk_no_metadata: None for target-group indexes -- not attributable to one group
    if nan_label is not None:
        stats["samples_dropped_missing_label"] = int(nan_label)
    if no_assembly is not None:
        stats["samples_dropped_missing_assembly"] = int(no_assembly)
    if on_disk_no_metadata is not None:
        stats["assemblies_excluded_missing_metadata"] = int(on_disk_no_metadata)
    stats["total_assemblies_written"] = int(len(metadata))

    stats_path = index_dir / f"{tag}_stats.json"
    with open(stats_path, "w") as f:
        json.dump(stats, f, indent=2)

    print(f"Wrote {index_type}: {len(metadata)} assemblies", file=sys.stderr)


def main():
    args = parse_args()

    metadata = pd.read_csv(args.metadata, low_memory=False)

    missing_cols = [c for c in [args.sample_col, args.group_label] if c not in metadata.columns]
    if missing_cols:
        sys.exit(f"Error: column(s) not found in metadata: {', '.join(missing_cols)}")

    nan_label = metadata[args.group_label].isna().sum()
    metadata = metadata.dropna(subset=[args.group_label])

    # sort by label lexicographically, then by sample ID within each group
    metadata = metadata.sort_values([args.group_label, args.sample_col], key=lambda col: col.astype(str))

    # Match assemblies to metadata
    assembly_input = args.assembly_dir if args.assembly_dir else args.assembly_paths
    metadata, no_assembly_by_group, on_disk_no_metadata = match_assemblies(
        metadata, assembly_input, args.assembly_suffix, args.sample_col, args.group_label
    )

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Count drops for stats
    no_assembly_count = sum(no_assembly_by_group.values())
    on_disk_no_metadata_count = len(on_disk_no_metadata)

    # Always write species-wide index
    write_index_files(
        metadata,
        output_dir,
        "species",
        None,
        args.group_label,
        args.sample_col,
        nan_label,
        no_assembly_count,
        on_disk_no_metadata_count,
    )

    # If target-groups provided, validate and write separate index per group
    target_groups = [g.strip() for g in args.target_groups.split(",") if g.strip()]
    if target_groups:
        found_groups = set(metadata[args.group_label].astype(str).unique())
        missing = set(target_groups) - found_groups
        if missing:
            print(f"Warning: target groups not found in metadata: {missing}", file=sys.stderr)

        for group in target_groups:
            group_metadata = metadata[metadata[args.group_label].astype(str) == group]
            if len(group_metadata) == 0:
                continue
            # this group's own no_assembly count (nan_label/on_disk_no_metadata don't apply per group)
            write_index_files(
                group_metadata,
                output_dir,
                "target_group",
                group,
                args.group_label,
                args.sample_col,
                no_assembly=no_assembly_by_group.get(group, 0),
            )


if __name__ == "__main__":
    main()
