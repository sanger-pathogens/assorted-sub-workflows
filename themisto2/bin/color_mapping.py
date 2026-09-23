#!/usr/bin/env python3
"""
Build Themisto2's ``--file-colors`` input for a species-wide colour index.

Reads one metadata table (one row per genome) and writes, into --output_dir,
three files prefixed with the species name (--species-name):

  <species>_file_colors_input.txt  assembly paths, one per line, grouped by label
  <species>_label_mapping.tsv      Sample_ID -> label, in index (colour-ID) order
  <species>_stats.json             summary counts, incl. assemblies per group

Every cell is read as text, so labels keep exactly the text in the file ("3"
never becomes "3.0"). Headers and the sample/label columns are stripped of
surrounding whitespace. A label that is empty or matches a missing value
(--label-missing, case-insensitive) becomes "unclassified".

Each label then goes through these rules, in order; the first that applies wins:

  1. --label-map    exact match on the raw label -> that map's group (final).
  2. missing        empty, or in --label-missing -> "unclassified".
  3. --label-multi  labels containing ";" (e.g. GPS merge history "1215;5"):
                    keep (as written) | smallest (lowest whole number, "5") |
                    unclassified.

Every change is listed in stats.json (label_changes) and summarised on stderr.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
from pathlib import Path

import pandas as pd

UNCLASSIFIED = "unclassified"
MULTI_SEP = ";"
MULTI_POLICIES = ("keep", "smallest", "unclassified")

# Group-label values treated as missing (-> UNCLASSIFIED), matched case-insensitively
# after stripping whitespace. --label-missing replaces this list; an empty label is
# always missing.
DEFAULT_MISSING = (
    "NA",
    "N/A",
    "#N/A",
    "NaN",
    "null",
    "none",
    "unknown",
    "missing",
    "-",
    "?",
    ".",
    "not applicable",
    "not available",
    "not collected",
    "not provided",
)

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
    p.add_argument(
        "--label-missing",
        default=None,
        help="'|'-separated group-label values to treat as missing (-> unclassified), matched "
        "case-insensitively, e.g. 'NA|unknown|not applicable'. Replaces the default list: "
        + ", ".join(repr(v) for v in DEFAULT_MISSING)
        + ". An empty label is always missing.",
    )
    p.add_argument(
        "--label-map",
        default=None,
        help="TSV with columns raw_label and group. An exact match on the raw label sets the final "
        "group, and no other rule touches it. Use group 'unclassified' to send a label to background.",
    )
    p.add_argument(
        "--label-multi",
        choices=MULTI_POLICIES,
        default="keep",
        help="What to do with labels containing ';': keep them as written (default), resolve to the "
        "smallest whole number ('1215;5' -> '5'), or send them to unclassified.",
    )
    p.add_argument("--output_dir", required=True, help="Directory to write the output files into.")
    return p.parse_args()


def parse_missing(raw: str | None) -> tuple[list[str], str]:
    """(missing values, 'built-in list' | '--label-missing') from --label-missing."""
    if raw is None:
        return list(DEFAULT_MISSING), "built-in list"
    return [v.strip() for v in raw.split("|") if v.strip()], "--label-missing"


def load_label_map(path: str, missing_set: set[str]) -> dict[str, str]:
    """raw_label -> group from a two-column TSV; exits listing every problem found."""
    label_map: dict[str, str] = {}
    errors: list[str] = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        reader.fieldnames = [(f or "").strip() for f in (reader.fieldnames or [])]
        if not {"raw_label", "group"} <= set(reader.fieldnames):
            sys.exit(f"Error: --label-map {path} needs columns 'raw_label' and 'group' (got {reader.fieldnames})")
        for line_no, row in enumerate(reader, start=2):
            raw, group = (row["raw_label"] or "").strip(), (row["group"] or "").strip()
            if not raw or not group:
                errors.append(f"line {line_no}: blank raw_label or group")
            elif group.casefold() in missing_set:
                errors.append(
                    f"line {line_no}: group '{group}' is a missing value -- use '{UNCLASSIFIED}' to send "
                    f"'{raw}' to background"
                )
            elif raw in label_map:
                errors.append(f"line {line_no}: raw_label '{raw}' is listed more than once")
            else:
                label_map[raw] = group
    if errors:
        sys.exit(f"Error: invalid --label-map {path}:\n  " + "\n  ".join(errors))
    return label_map


def resolve_label(raw: str, label_map: dict, missing_set: set[str], multi: str) -> tuple[str, str | None]:
    """(final label, rule that changed it: 'label_map' | 'label_missing' | 'label_multi' | None).

    Raises ValueError for label_multi=smallest on a label whose parts aren't all
    whole numbers -- never guess.
    """
    if raw in label_map:
        return label_map[raw], "label_map"
    if raw == "" or raw.casefold() in missing_set:
        return UNCLASSIFIED, "label_missing"
    if MULTI_SEP in raw and multi != "keep":
        if multi == "unclassified":
            return UNCLASSIFIED, "label_multi"
        parts = [p.strip() for p in raw.split(MULTI_SEP)]
        if not all(p.isascii() and p.isdigit() for p in parts):
            raise ValueError(raw)
        return str(min(int(p) for p in parts)), "label_multi"
    return raw, None


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

    sample_col, group_label = args.sample_col.strip(), args.group_label.strip()
    missing_values, missing_source = parse_missing(args.label_missing)
    missing_set = {v.casefold() for v in missing_values}

    # All text, no pandas NA guessing: labels are exactly what's in the file.
    metadata = pd.read_csv(args.metadata, dtype=str, keep_default_na=False)
    metadata.columns = [str(c).strip() for c in metadata.columns]

    missing = [c for c in (sample_col, group_label) if c not in metadata.columns]
    if missing:
        sys.exit(f"Error: column(s) not found in metadata: {', '.join(missing)}")

    metadata[sample_col] = metadata[sample_col].str.strip()
    metadata[group_label] = metadata[group_label].str.strip()
    duplicates = int(metadata[sample_col].duplicated().sum())
    if duplicates:
        metadata = metadata.drop_duplicates(subset=sample_col)
        print(f"warning: ignored {duplicates} row(s) with a duplicate {sample_col}", file=sys.stderr)

    label_map = load_label_map(args.label_map, missing_set) if args.label_map else {}

    # Resolve each distinct raw label once. Genomes whose label is missing aren't
    # dropped -- they're relabelled unclassified and kept.
    resolved, bad_multi = {}, []
    for raw in metadata[group_label].unique():
        try:
            resolved[raw] = resolve_label(raw, label_map, missing_set, args.label_multi)
        except ValueError:
            bad_multi.append(raw)
    if bad_multi:
        sys.exit(
            f"Error: --label-multi smallest needs every '{MULTI_SEP}'-separated part to be a whole number; "
            f"these labels aren't: {', '.join(repr(b) for b in sorted(bad_multi))}. "
            "Fix them in the metadata or list them in --label-map."
        )
    metadata["_raw_label"] = metadata[group_label]
    metadata["_rule"] = metadata["_raw_label"].map(lambda r: resolved[r][1])
    metadata[group_label] = metadata["_raw_label"].map(lambda r: resolved[r][0])
    blank_label = metadata[group_label].eq(UNCLASSIFIED) & metadata["_rule"].notna()
    label_map_unmatched = sorted(set(label_map) - set(resolved))

    # Match each metadata row to an assembly by normalised filename.
    assemblies, dead_paths = load_assemblies(args.assembly_dir or args.assembly_paths)
    if dead_paths:
        print(
            f"warning: {len(dead_paths)} assembly path(s) don't point at a file, e.g. "
            + ", ".join(dead_paths[:3]),
            file=sys.stderr,
        )
    by_safe_name = index_by_safe_name(assemblies)
    metadata["_want"] = (metadata[sample_col] + args.assembly_suffix).map(fs_safe)
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
            sample_col: [strip_suffix(Path(p).name, args.assembly_suffix) for p in orphan_paths],
            group_label: UNCLASSIFIED,
            "file_path": orphan_paths,
        }
    )

    written = pd.concat(
        [kept[[sample_col, group_label, "file_path"]], orphan_rows], ignore_index=True
    )
    # sort by label, then sample ID within each label -- fixes colour-ID order
    written = written.sort_values(
        [group_label, sample_col], key=lambda c: c.astype(str)
    )

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    written["file_path"].to_csv(out / f"{prefix}_file_colors_input.txt", index=False, header=False)
    (
        written[[sample_col, group_label]]
        .rename(columns={sample_col: "Sample_ID", group_label: "label"})
        .to_csv(out / f"{prefix}_label_mapping.tsv", index=False, sep="\t")
    )

    per_group = {
        str(label): int(n)
        for label, n in written[group_label].value_counts().sort_index().items()
    }
    rewritten = kept[kept[group_label] != kept["_raw_label"]]
    label_changes = [
        {"raw_label": raw, "new_label": final, "changed_by": rule, "genomes": int(n)}
        for (raw, final, rule), n in rewritten.groupby(["_raw_label", group_label, "_rule"]).size().items()
    ]
    label_changes.sort(key=lambda r: (-r["genomes"], r["raw_label"]))

    stats = {
        "species": prefix,
        "group_label_column": group_label,
        "assemblies_total": len(written) + n_dropped,
        "assemblies_written": len(written),
        "assemblies_dropped_fasta_not_found": n_dropped,
        "assembly_paths_missing_file": len(dead_paths),
        "relabelled_unclassified": n_relabelled,
        "assemblies_without_metadata_row": len(orphan_paths),
        "assemblies_per_group": per_group,
        "label_settings": {
            "label_missing": missing_values,
            "label_missing_source": missing_source,
            "label_multi": args.label_multi,
            "label_map": args.label_map,
        },
        "label_changes": label_changes,
        "label_map_unmatched": label_map_unmatched,
        "note": (
            "Genome counts first, then how group labels were cleaned. 'unclassified' genomes are "
            "metadata rows whose label was blank or a missing value (or was sent there by a label "
            "rule), plus assembly files with no metadata row. "
            "They stay in the index as one group. "
            "assemblies_dropped_fasta_not_found are metadata rows whose assembly FASTA wasn't found."
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
        summary += f"; {n_relabelled} missing or rewritten label -> {UNCLASSIFIED}"
    print(summary, file=sys.stderr)
    if label_changes:
        print(
            f"{prefix}: {len(label_changes)} label(s) changed (see label_changes in stats.json):",
            file=sys.stderr,
        )
        for r in label_changes[:20]:
            print(
                f"  '{r['raw_label']}' -> '{r['new_label']}' (changed by {r['changed_by']}, {r['genomes']} genomes)",
                file=sys.stderr,
            )
        if len(label_changes) > 20:
            print(f"  ... and {len(label_changes) - 20} more", file=sys.stderr)
    if label_map_unmatched:
        print(
            f"warning: {len(label_map_unmatched)} --label-map raw_label(s) match no metadata label: "
            + ", ".join(repr(u) for u in label_map_unmatched[:10]),
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
