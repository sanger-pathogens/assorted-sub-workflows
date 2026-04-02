#!/usr/bin/env python3

import argparse
import logging
import re
from pathlib import Path

import pandas as pd

TAXONOMIC_RANK = {
    "d": "domain",
    "k": "kingdom",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
    "t": "taxon"
}


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sylph_prof_report", help="Report file (TSV) with `sylph profile` output format.", type=Path, required=True)
    parser.add_argument("--sylphtax_report", help="One or more report files (TSV) with sylph-tax taxprof output format.", type=Path, nargs="+", required=True)
    parser.add_argument("--outdir", type=Path, default=".")
    parser.add_argument("--taxonomic_group", help="Taxonomic rank by which to group report", type=str, choices=list(TAXONOMIC_RANK.values()), default="species")
    parser.add_argument("--prefix", help="Prefix for output report filenames", type=str)
    return parser.parse_args()


def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

def load_data(report: Path) -> pd.DataFrame:
    df = pd.read_csv(report, sep="\t", comment="#")
    return df


def load_taxonomy_data(reports) -> pd.DataFrame:
    frames = [load_data(report) for report in reports]
    return pd.concat(frames, ignore_index=True)

def extract_genbank_id(df: pd.DataFrame, column: str) -> pd.DataFrame:
    extract_pattern = {
        "clade_name": r"\|t__(.*)",
        "Genome_file": r".*/(.*)_genomic.fna.gz"
    }
    new_df = df.copy()
    try:
        pattern = extract_pattern[column]
    except KeyError:
        raise ValueError(f"Given column '{column}' is not supported.")
    new_df["genome_id"] = new_df[column].str.extract(pattern)
    return new_df

def join_taxonomy_to_profile(sylph_prof_df: pd.DataFrame, sylphtax_df: pd.DataFrame, join_column: str) -> pd.DataFrame:
    return sylph_prof_df.merge(sylphtax_df, on=join_column, how="left")

def split_taxonomy(df: pd.DataFrame, taxonomy_column: str) -> pd.DataFrame:
    new_df = df.copy()
    for label, rank in TAXONOMIC_RANK.items():
        pattern = fr"{label}__([^|]*)"
        new_df[rank] = new_df[taxonomy_column].str.extract(pattern)
    return new_df

def normalize_taxon_name(taxon: str) -> str:
    return re.sub(r"\s+", "_", taxon.lower())

def main():
    args = parse_args()
    setup_logging()

    outdirs = [args.outdir / "reports", args.outdir / "refs"]
    for outdir in outdirs:
        outdir.mkdir(parents=True, exist_ok=True)

    sylph_prof_df = extract_genbank_id(load_data(args.sylph_prof_report), column="Genome_file")
    sylphtax_df = extract_genbank_id(load_taxonomy_data(args.sylphtax_report), column="clade_name")
    sylphtax_df = sylphtax_df.drop_duplicates(subset=["genome_id", "clade_name"])

    columns_to_drop = ["relative_abundance", "sequence_abundance", "ANI (if strain-level)", "Coverage (if strain-level)"]
    merged_df = sylph_prof_df.merge(sylphtax_df, on="genome_id", how="left").drop(columns_to_drop, axis=1)

    merged_df = split_taxonomy(merged_df, taxonomy_column="clade_name")

    dfs = {taxon_name: df for taxon_name, df in merged_df.groupby(args.taxonomic_group)}

    for taxon_name, df in dfs.items():
        report_prefix = normalize_taxon_name(taxon_name)
        if args.prefix:
            report_prefix = f"{args.prefix}_{report_prefix}"
        df.to_csv(args.outdir / "reports" / f"{report_prefix}.tsv", sep="\t", index=False)
        # Ensure we have a unique list of genomes as output
        refs = pd.Series(df["Genome_file"].unique())
        refs.to_csv(args.outdir / "refs" / f"{report_prefix}.txt", sep="\t", index=False, header=False)

if __name__ == "__main__":
    main()