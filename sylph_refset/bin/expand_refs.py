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
    parser.add_argument("--sylphtax_report", help="One or more report files (TSV) with sylph-tax taxprof output format.", type=Path, nargs="+", required=True)
    parser.add_argument("--taxonomy_data", help="Two column (TSV) file linking genome IDs to metaphlan-format taxonomy.", type=Path, required=True)
    parser.add_argument("--genome_to_file", help="Two column (TSV) file linking genome IDs to genome fasta filepaths.", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, default=".")
    parser.add_argument("--taxonomic_group", help="Taxonomic rank by which to group report.", type=str, choices=list(TAXONOMIC_RANK.values()), default="species")
    parser.add_argument("--remove_pattern", help="Pattern (regex) to remove some string from labels of the given taxonomic groups.", type=str)
    parser.add_argument("--prefix", help="Prefix for output report filenames", type=str)
    return parser.parse_args()


def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def load_data(report: Path, header=None) -> pd.DataFrame:
    df = pd.read_csv(report, sep="\t", comment="#", header=header)
    return df


def load_taxonomy_data(reports) -> pd.DataFrame:
    frames = [load_data(report, header=0) for report in reports]
    return pd.concat(frames, ignore_index=True)


def split_taxonomy(df: pd.DataFrame, taxonomy_column: str) -> pd.DataFrame:
    new_df = df.copy()
    for label, rank in TAXONOMIC_RANK.items():
        pattern = fr"{label}__([^|;]*)"
        new_df[rank] = new_df[taxonomy_column].str.extract(pattern)
    return new_df

def normalize_taxon_name(taxon: str) -> str:
    return re.sub(r"\s+", "_", taxon.lower())

def main():
    args = parse_args()
    setup_logging()

    args.outdir.mkdir(parents=True, exist_ok=True)

    sylphtax_df = load_taxonomy_data(args.sylphtax_report)
    sylphtax_df = split_taxonomy(sylphtax_df, taxonomy_column="clade_name")
    taxa = sylphtax_df[args.taxonomic_group].unique()
    # print(taxa)
    
    taxonomy_df = load_data(args.taxonomy_data)
    taxonomy_df = split_taxonomy(taxonomy_df, taxonomy_column=1)
    taxonomy_df_subset = taxonomy_df.loc[taxonomy_df[args.taxonomic_group].isin(taxa)]
    if args.remove_pattern:
        remove_pattern = re.compile(args.remove_pattern)
        taxonomy_df_subset[args.taxonomic_group] = taxonomy_df_subset[args.taxonomic_group].str.replace(remove_pattern, "", regex=True)
    # print(taxonomy_df_subset.head)
    taxonomy_specific_dfs = {taxon_name: df for taxon_name, df in taxonomy_df_subset.groupby(args.taxonomic_group)}

    genome_id_to_file = load_data(args.genome_to_file)
    
    for taxon_name, df in taxonomy_specific_dfs.items():
        report_prefix = normalize_taxon_name(taxon_name)
        genome_ids = df[0]
        genome_files = genome_id_to_file.loc[genome_id_to_file[0].isin(genome_ids), 1]

        if args.prefix:
            report_prefix = f"{args.prefix}_{report_prefix}"
        genome_files.to_csv(args.outdir / f"{report_prefix}.txt", sep="\t", index=False, header=False)


if __name__ == "__main__":
    main()
