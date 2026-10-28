#!/usr/bin/env python3
"""Concatenate per-batch qc_mags reports into one table with a sample_id column.

sample_id comes from the manifest (genome -> manifest ID map built before
batching), matched on preqc_genome_name.
"""

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

GENOME_COLUMN = "preqc_genome_name"
SAMPLE_COLUMN = "sample_id"
REPORT_SUFFIX = "_final_report.tsv"


def load_sample_map(path: Path) -> dict:
    """genome_name -> list of sample IDs (a name can repeat across samples)."""
    if path.stat().st_size == 0:
        return {}
    m = pd.read_csv(path, sep="\t", header=None, names=["sample_id", "genome_name"],
                    dtype=str, keep_default_na=False)
    return m.drop_duplicates().groupby("genome_name")["sample_id"].apply(list).to_dict()


def resolve(genome: str, batch_id: str, sample_map: dict):
    samples = sample_map.get(genome)
    if not samples:
        return None
    if len(samples) == 1:
        return samples[0]
    # Same genome file name in several samples: only distinguishable when each
    # sample is its own batch (batch_size 0), where the batch ID is the sample ID.
    if batch_id in samples:
        return batch_id
    return ";".join(sorted(samples))


def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reports", type=Path, nargs="+", required=True)
    ap.add_argument("--sample_map", type=Path, required=True,
                    help="TSV without header: sample_id <TAB> genome_name")
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    sample_map = load_sample_map(args.sample_map)
    frames = []
    for rp in sorted(args.reports):
        df = pd.read_csv(rp, sep="\t", dtype=str, keep_default_na=False)
        batch_id = rp.name.removesuffix(REPORT_SUFFIX)
        df.insert(1, SAMPLE_COLUMN, [resolve(g, batch_id, sample_map) for g in df[GENOME_COLUMN]])
        frames.append(df)
        logging.info(f"{rp.name}: {len(df)} genomes")

    combined = pd.concat(frames, ignore_index=True, sort=False)
    n_missing = int(combined[SAMPLE_COLUMN].isna().sum())
    if n_missing:
        logging.warning(f"{n_missing}/{len(combined)} genomes have no sample_id in the manifest map")
    n_ambig = int(combined[SAMPLE_COLUMN].fillna("").str.contains(";").sum())
    if n_ambig:
        logging.warning(f"{n_ambig} genome names occur in several samples; sample_id lists all of them")
    combined.to_csv(args.output, sep="\t", index=False)
    logging.info(f"Combined report: {len(combined)} genomes from {len(frames)} batch reports -> {args.output}")


if __name__ == "__main__":
    sys.exit(main())
