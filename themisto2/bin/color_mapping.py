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
surrounding whitespace. Each label then goes through these rules, in order:

  1. missing  empty, or one of MISSING_VALUES (case-insensitive) -> "unclassified".
  2. GPSC     only when --group-label is GPSC (any case): a label containing ";"
              (GPS merge history, e.g. "1215;5" or "GPSC1215;GPSC5") becomes its
              smallest number, keeping the label's GPSC prefix if it has one
              ("GPSC1215;5" -> "GPSC5"). The one exception is 235 with 9, in any
              order or prefix form, which is kept as its own group "235_9"
              ("GPSC235_9" with the prefix): it's a mixture of GPSC9 and GPSC235,
              but current evidence doesn't say to merge them. A ";" label with a
              part that isn't a whole number stops the run, listing every bad label.
              For any other --group-label, ";" labels are left as written.

Every change is listed in stats.json (label_changes) and summarised on stderr.

Genomes whose label is "unclassified" (a missing value, a label that reads
"unclassified", or an assembly with no metadata row) are always left out of the
index, so markers are never checked against them. They're listed, with the
reason, in <species>_dropped_unclassified.tsv. The run stops if nothing is left.
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
MULTI_SEP = ";"

# Group-label values treated as missing (-> UNCLASSIFIED), matched case-insensitively
# after stripping whitespace. An empty label is always missing.
MISSING_VALUES = (
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
MISSING_SET = {v.casefold() for v in MISSING_VALUES}

# One part of a GPSC ";" label: an optional GPSC prefix and a whole number.
GPSC_PART = re.compile(r"(?i)(gpsc)?(\d+)")
# GPSC235;9 is a mixture of GPSC9 and GPSC235, but there's no evidence to merge the
# two lineages, so it's kept as its own group instead of taking the smallest number.
GPSC_UNMERGED = {9, 235}

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


def resolve_label(raw: str, gpsc: bool) -> tuple[str, str | None]:
    """(final label, rule that changed it: 'missing_value' | 'gpsc_multi' | None).

    Raises ValueError for a GPSC ";" label whose parts aren't all whole numbers
    (optionally GPSC-prefixed) -- never guess.
    """
    if raw == "" or raw.casefold() in MISSING_SET:
        return UNCLASSIFIED, "missing_value"
    if not (gpsc and MULTI_SEP in raw):
        return raw, None
    parts = [GPSC_PART.fullmatch(p.strip()) for p in raw.split(MULTI_SEP)]
    if not all(parts):
        raise ValueError(raw)
    numbers = {int(m.group(2)) for m in parts}
    prefix = parts[0].group(1) or ""
    if numbers == GPSC_UNMERGED:
        return f"{prefix}235_9", "gpsc_multi"
    return f"{prefix}{min(numbers)}", "gpsc_multi"


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

    # Resolve each distinct raw label once.
    gpsc = group_label.casefold() == "gpsc"
    resolved, bad_multi = {}, []
    for raw in metadata[group_label].unique():
        try:
            resolved[raw] = resolve_label(raw, gpsc)
        except ValueError:
            bad_multi.append(raw)
    if bad_multi:
        sys.exit(
            f"Error: GPSC labels containing '{MULTI_SEP}' must be whole numbers (optionally GPSC-prefixed) "
            f"in every part; these aren't: {', '.join(repr(b) for b in sorted(bad_multi))}. "
            "Fix them in the metadata."
        )
    metadata["_raw_label"] = metadata[group_label]
    metadata["_rule"] = metadata["_raw_label"].map(lambda r: resolved[r][1])
    metadata[group_label] = metadata["_raw_label"].map(lambda r: resolved[r][0])
    # A label that already reads "unclassified" (any case) is unclassified too.
    labelled_uncl = metadata[group_label].str.casefold().eq(UNCLASSIFIED)
    metadata.loc[labelled_uncl, group_label] = UNCLASSIFIED
    blank_label = metadata["_rule"].eq("missing_value")

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

    # Assembly files no metadata row claimed are "unclassified" (and dropped below).
    claimed = set(kept["_want"])
    orphan_paths = sorted(p for safe, p in by_safe_name.items() if safe not in claimed)
    orphan_rows = pd.DataFrame(
        {
            sample_col: [strip_suffix(Path(p).name, args.assembly_suffix) for p in orphan_paths],
            group_label: UNCLASSIFIED,
            "file_path": orphan_paths,
            "_raw_label": "",
            "_reason": "no_metadata_row",
        }
    )
    # Why each genome is unclassified, for the dropped list: a missing value, or a
    # metadata label that literally reads "unclassified".
    kept["_reason"] = kept["_rule"].where(kept["_rule"].eq("missing_value"), "labelled_unclassified")

    cols = [sample_col, group_label, "file_path", "_raw_label", "_reason"]
    written = pd.concat([kept[cols], orphan_rows[cols]], ignore_index=True)

    is_uncl = written[group_label] == UNCLASSIFIED
    dropped_unclassified, written = written[is_uncl], written[~is_uncl]
    if written.empty:
        sys.exit(
            f"Error: no genomes left for {prefix} -- all {len(dropped_unclassified)} are "
            f"{UNCLASSIFIED}. Check --group-label and the metadata labels."
        )
    # sort by label, then sample ID within each label -- fixes colour-ID order
    written = written.sort_values(
        [group_label, sample_col], key=lambda c: c.astype(str)
    )

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    dropped_path = out / f"{prefix}_dropped_unclassified.tsv"
    (
        dropped_unclassified[[sample_col, "_raw_label", "_reason"]]
        .rename(columns={sample_col: "Sample_ID", "_raw_label": "raw_label", "_reason": "reason"})
        .sort_values(["reason", "Sample_ID"])
        .to_csv(dropped_path, index=False, sep="\t")
    )

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
    rewritten = kept[kept["_rule"].notna()]
    label_changes = [
        {"raw_label": raw, "new_label": final, "changed_by": rule, "genomes": int(n)}
        for (raw, final, rule), n in rewritten.groupby(["_raw_label", group_label, "_rule"]).size().items()
    ]
    label_changes.sort(key=lambda r: (-r["genomes"], r["raw_label"]))

    stats = {
        "species": prefix,
        "group_label_column": group_label,
        "assemblies_total": len(written) + n_dropped + len(dropped_unclassified),
        "assemblies_written": len(written),
        "assemblies_dropped_fasta_not_found": n_dropped,
        "assemblies_dropped_unclassified": len(dropped_unclassified),
        "assembly_paths_missing_file": len(dead_paths),
        "relabelled_unclassified": n_relabelled,
        "assemblies_without_metadata_row": len(orphan_paths),
        "assemblies_per_group": per_group,
        "missing_values": list(MISSING_VALUES),
        "label_changes": label_changes,
        "note": (
            "Genome counts first, then how group labels were cleaned. 'unclassified' genomes are "
            "metadata rows whose label was blank, a missing value or read 'unclassified', plus "
            "assembly files with no metadata row. They're always left out of the index: counted in "
            "assemblies_dropped_unclassified and listed in <species>_dropped_unclassified.tsv. "
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
    if len(dropped_unclassified):
        summary += f"; {len(dropped_unclassified)} {UNCLASSIFIED} dropped from the index"
    if n_relabelled:
        summary += f"; {n_relabelled} missing label(s) -> {UNCLASSIFIED}"
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
    if len(dropped_unclassified):
        print(
            f"warning: {len(dropped_unclassified)} {UNCLASSIFIED} genome(s) left out of the index; "
            f"markers are not checked against them. See {dropped_path.name}",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
