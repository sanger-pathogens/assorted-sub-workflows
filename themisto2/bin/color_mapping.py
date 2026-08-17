#!/usr/bin/env python3
"""
Generic color mapping script for Themisto2 --file-colors input.

Reads metadata CSV, sorts by label column, matches samples to assemblies,
and prepares index input files organized in subdirectories:
  - index_species/ : whole-species index files
  - index_target_group/<group>/ : per-group indexes (if --target-groups provided)
"""

import argparse
import json
import os
import sys
from pathlib import Path

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build Themisto2 file-colors input from a metadata table."
    )
    parser.add_argument(
        "--metadata", required=True,
        help="Path to metadata CSV file."
    )
    parser.add_argument(
        "--sample-col", default="Sample_ID",
        help="Column name in metadata containing sample identifiers (default: Sample_ID)."
    )
    parser.add_argument(
        "--label-col", required=True,
        help="Column name to sort samples by (lexicographically) and write to label_mapping.tsv."
    )
    parser.add_argument(
        "--target-groups", default="",
        help="Comma-separated label(s) to build separate lineage-scoped indexes for, e.g. GPSC1,GPSC2. Leave empty for species-wide only."
    )
    asm = parser.add_mutually_exclusive_group(required=True)
    asm.add_argument(
        "--assembly-dir",
        help="Directory containing assembly FASTA files."
    )
    asm.add_argument(
        "--assembly-paths",
        help="Text file listing one assembly path per line."
    )
    parser.add_argument(
        "--assembly-suffix", default=".contigs.fasta",
        help="Suffix appended to Sample_ID to form the assembly filename (default: .contigs.fasta). Only used with --assembly-dir."
    )
    parser.add_argument(
        "--output_dir", required=True,
        help="Output directory for index subdirectories."
    )
    return parser.parse_args()


def match_assemblies(metadata, assembly_input, assembly_suffix, sample_col):
    """
    Match metadata samples to assembly files on disk.
    
    Args:
        metadata: DataFrame with sample identifiers
        assembly_input: directory path or file listing assembly paths
        assembly_suffix: suffix to append to sample IDs
        sample_col: column name containing sample identifiers
    
    Returns:
        metadata: updated with file_path column
        no_assembly: boolean mask of unmatched samples
        on_disk_no_metadata: set of assemblies with no metadata entry
    """
    assembly_input_path = Path(assembly_input)
    
    if assembly_input_path.is_dir():
        # Resolve to canonical absolute path so file_colors_input.txt can be read
        # by GGCAT in a different task directory (symlinks won't resolve there)
        assembly_dir = assembly_input_path.resolve()
        assembly_files = set(os.listdir(assembly_dir))
        metadata["filename"] = metadata[sample_col].astype(str) + assembly_suffix
        no_assembly = ~metadata["filename"].isin(assembly_files)
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
        on_disk_no_metadata = set(path_lookup.keys()) - set(metadata["filename"])
        metadata = metadata[metadata["filename"].isin(path_lookup)].reset_index(drop=True)
        metadata["file_path"] = metadata["filename"].apply(lambda f: path_lookup[f])
    
    return metadata, no_assembly, on_disk_no_metadata


def write_index_files(metadata, output_dir, index_type, group_name, label_col, sample_col, nan_label, no_assembly, on_disk_no_metadata):
    """
    Write index files to appropriate subdirectory.
    
    Args:
        metadata: filtered DataFrame
        output_dir: Path object for output directory
        index_type: "species" or "target_group"
        group_name: group label (e.g., "GPSC1"), None for species-wide
        label_col: column name for labels
        sample_col: column name for samples
        nan_label: count of samples dropped for no label
        no_assembly: count of samples dropped for no assembly on disk
        on_disk_no_metadata: count of assemblies on disk with no metadata entry
    """
    if index_type == "species":
        index_dir = output_dir / "index_species"
    else:
        index_dir = output_dir / "index_target_group" / group_name
    
    index_dir.mkdir(parents=True, exist_ok=True)
    
    metadata["file_path"].to_csv(
        index_dir / "file_colors_input.txt", index=False, header=False
    )

    label_mapping = (
        metadata[[sample_col, label_col]]
        .drop_duplicates()
        .rename(columns={sample_col: "Sample_ID", label_col: "label"})
    )
    label_mapping.to_csv(index_dir / "label_mapping.tsv", index=False, sep="\t")

    stats = {
        "index_type": index_type,
        "target_group": group_name or "species_wide",
        "metadata_column": label_col,
        "samples_no_label_dropped": int(nan_label),
        "samples_no_assembly_dropped": int(no_assembly),
        "assemblies_no_metadata_excluded": int(on_disk_no_metadata),
        "total_assemblies_written": int(len(metadata)),
    }
    
    stats_path = index_dir / "stats.json"
    with open(stats_path, "w") as f:
        json.dump(stats, f, indent=2)

    print(f"Wrote {index_type}: {len(metadata)} assemblies", file=sys.stderr)


def main():
    args = parse_args()

    metadata = pd.read_csv(args.metadata, low_memory=False)

    missing_cols = [c for c in [args.sample_col, args.label_col] if c not in metadata.columns]
    if missing_cols:
        sys.exit(f"Error: column(s) not found in metadata: {', '.join(missing_cols)}")

    nan_label = metadata[args.label_col].isna().sum()
    metadata = metadata.dropna(subset=[args.label_col])

    # sort by label lexicographically, then by sample ID within each group
    metadata = metadata.sort_values(
        [args.label_col, args.sample_col],
        key=lambda col: col.astype(str)
    )

    # Match assemblies to metadata
    assembly_input = args.assembly_dir if args.assembly_dir else args.assembly_paths
    metadata, no_assembly, on_disk_no_metadata = match_assemblies(
        metadata, assembly_input, args.assembly_suffix, args.sample_col
    )

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Count drops for stats
    no_assembly_count = no_assembly.sum()
    on_disk_no_metadata_count = len(on_disk_no_metadata)

    # Always write species-wide index
    write_index_files(metadata, output_dir, "species", None, args.label_col, args.sample_col, 
                nan_label, no_assembly_count, on_disk_no_metadata_count)

    # If target-groups provided, validate and write separate index per group
    target_groups = [g.strip() for g in args.target_groups.split(",") if g.strip()]
    if target_groups:
        found_groups = set(metadata[args.label_col].astype(str).unique())
        missing = set(target_groups) - found_groups
        if missing:
            print(f"Warning: target groups not found in metadata: {missing}", file=sys.stderr)
        
        for group in target_groups:
            group_metadata = metadata[metadata[args.label_col].astype(str) == group]
            if len(group_metadata) == 0:
                continue
            write_index_files(group_metadata, output_dir, "target_group", group, args.label_col, args.sample_col,
                        0, 0, 0)


if __name__ == "__main__":
    main()