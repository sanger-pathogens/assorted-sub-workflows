#!/usr/bin/env python3

"""
Clean and check input genome FASTA files before QC.

For each input file (plain or gzipped), scans for problems and appends any fixes/warnings to a TSV:
    - removes GenBank-style '//' record terminator lines
    - removes blank (whitespace-only) lines
    - warns (does not fix) on sequence lines containing non-IUPAC nucleotide characters

Files with lines to remove get a cleaned copy (same name and compression) in the output
directory; all other files are symlinked there unchanged, to avoid copying every genome.

Fixes TSV columns (no header): batch, file name, fix type, line count

Usage:
    check_input_fasta.py \
        --batch <batch_id> \
        --fixes <output/fixes.tsv> \
        --outdir <checked/> \
        <fasta> [<fasta> ...]
"""

import argparse
import gzip
import re
from functools import partial
from pathlib import Path

GZIP_MAGIC = b"\x1f\x8b"
WHITESPACE = " \t\r\v\f"
# Any character that is not an IUPAC nucleotide code
NON_NUCLEOTIDE = re.compile(r"[^ACGTNRYKMSWBDHVacgtnrykmswbdhv]")


def main():
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    with open(args.fixes, "w") as fixes:
        for fasta in sorted(args.fastas):
            slash, blank, bad = check_fasta(fasta, args.outdir / fasta.name)
            for fix, count in (
                ("removed_double_slash_lines", slash),
                ("removed_blank_lines", blank),
                ("WARNING_non_nucleotide_characters", bad),
            ):
                if count:
                    fixes.write(f"{args.batch}\t{fasta.name}\t{fix}\t{count}\n")


def parse_args():
    parser = argparse.ArgumentParser(description="Clean and check input genome FASTA files.")
    parser.add_argument("--batch", required=True, help="Batch ID written to the fixes TSV")
    parser.add_argument("--fixes", required=True, type=Path, help="Output fixes/warnings TSV")
    parser.add_argument("--outdir", required=True, type=Path, help="Directory for cleaned FASTAs")
    parser.add_argument("fastas", nargs="+", type=Path, help="Input FASTA files (optionally gzipped)")
    return parser.parse_args()


def is_gzipped(path):
    with open(path, "rb") as f:
        return f.read(2) == GZIP_MAGIC


def check_fasta(src, dest):
    """Scan src; write a cleaned copy to dest if lines need removing, else symlink it.
    Returns (slash, blank, bad) line counts."""
    opener = partial(gzip.open, compresslevel=6) if is_gzipped(src) else open
    text_mode = dict(encoding="latin-1", newline="")

    with opener(src, "rt", **text_mode) as fin:
        slash, blank, bad = clean_lines(fin)

    if slash or blank:
        with opener(src, "rt", **text_mode) as fin, opener(dest, "wt", **text_mode) as fout:
            clean_lines(fin, fout)
    else:
        dest.symlink_to(src.resolve())
    return slash, blank, bad


def clean_lines(fin, fout=None):
    """Count (and, if fout is given, write out) lines with '//' and blank lines dropped.
    Returns (slash, blank, bad) line counts."""
    slash = blank = bad = 0
    for line in fin:
        line = line.rstrip("\n")
        if line == "//":
            slash += 1
            continue
        if not line.strip(WHITESPACE):
            blank += 1
            continue
        if not line.startswith(">") and NON_NUCLEOTIDE.search(line):
            bad += 1
        if fout is not None:
            fout.write(line + "\n")
    return slash, blank, bad


if __name__ == "__main__":
    main()
