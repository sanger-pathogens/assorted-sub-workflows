#!/usr/bin/env python3
"""
Generic color mapping script for Themisto2 --file-colors input.

Reads metadata CSV, sorts by group label column, matches samples to
assemblies, and writes the species-wide index input files under
index_species/ (file_colors_input.txt, label_mapping.tsv, stats.json).

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
        help="Column name to sort samples by (lexicographically) and write to label_mapping.tsv.",
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
    group_label: str,
    sample_col: str,
    nan_label: int,
    no_assembly: int,
    on_disk_no_metadata: int,
) -> None:
    """Write index_species/{species_file_colors_input.txt, species_label_mapping.tsv, species_stats.json}.

    Args:
        metadata: All kept rows (already filtered/sorted), with a ``file_path``
            column plus ``sample_col`` and ``group_label``.
        output_dir: Root output dir; files go under ``index_species/``.
        group_label: Metadata column used for labelling (label_mapping.tsv's ``label``).
        sample_col: Metadata column holding sample identifiers.
        nan_label: Count of samples dropped for a missing group label.
        no_assembly: Count of samples dropped for having no assembly on disk.
        on_disk_no_metadata: Count of assemblies on disk with no metadata row.
    """
    index_dir = output_dir / "index_species"
    index_dir.mkdir(parents=True, exist_ok=True)

    metadata["file_path"].to_csv(index_dir / "species_file_colors_input.txt", index=False, header=False)

    label_mapping = (
        metadata[[sample_col, group_label]]
        .drop_duplicates()
        .rename(columns={sample_col: "Sample_ID", group_label: "label"})
    )
    label_mapping.to_csv(index_dir / "species_label_mapping.tsv", index=False, sep="\t")

    stats = {
        "index_type": "species",
        "metadata_column": group_label,
        "samples_dropped_missing_label": int(nan_label),
        "samples_dropped_missing_assembly": int(no_assembly),
        "assemblies_excluded_missing_metadata": int(on_disk_no_metadata),
        "total_assemblies_written": int(len(metadata)),
    }
    with open(index_dir / "species_stats.json", "w") as f:
        json.dump(stats, f, indent=2)

    print(f"Wrote species index: {len(metadata)} assemblies", file=sys.stderr)


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

    write_index_files(
        metadata,
        output_dir,
        args.group_label,
        args.sample_col,
        nan_label,
        no_assembly_count,
        on_disk_no_metadata_count,
    )


if __name__ == "__main__":
    main()
