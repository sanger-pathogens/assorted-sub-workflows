#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd

SUMMARY_COLUMNS = ["species_name", "removed_reference_count", "removed_genome_ids"]


def parse_args():
    p = argparse.ArgumentParser(description="Filter a combined sylph TSV report using ANI and effective coverage thresholds.")
    p.add_argument("--input", type=Path, required=True, help="Combined sylph report TSV.")
    p.add_argument(
        "--taxonomy-data",
        type=Path,
        help="Two-column taxonomy TSV mapping genome IDs to GTDB taxonomy strings.",
    )
    p.add_argument("--ani", type=float, required=True, help="Minimum ANI threshold.")
    p.add_argument("--cov", type=float, required=True, help="Minimum Eff_cov threshold.")
    p.add_argument("--ani-column", default="Adjusted_ANI")
    p.add_argument("--cov-column", default="Eff_cov")
    p.add_argument("--out-report", type=Path, required=True, help="Filtered TSV output path.")
    p.add_argument("--out-summary", type=Path, help="Optional summary TSV output path.")
    p.add_argument(
        "--out-references",
        type=Path,
        help="Optional file containing unique Genome_file entries from the filtered report.",
    )
    return p.parse_args()


def load_species_by_genome_id(path):
    if not path:
        return {}

    try:
        taxonomy_df = pd.read_csv(
            path,
            sep="\t",
            header=None,
            usecols=[0, 1],
            names=["genome_id", "taxonomy"],
        )
    except pd.errors.EmptyDataError:
        return {}

    taxonomy_df["species_name"] = (
        taxonomy_df["taxonomy"]
        .fillna("")
        .str.extract(r"s__([^|;]*)", expand=False)
        .fillna("unknown_species")
        .replace("", "unknown_species")
    )
    return dict(zip(taxonomy_df["genome_id"], taxonomy_df["species_name"]))


def main():
    args = parse_args()
    df = pd.read_csv(args.input, sep="\t")

    missing_columns = sorted({"Genome_file", args.ani_column, args.cov_column}.difference(df.columns))
    if missing_columns:
        raise ValueError(f"Input report is missing required columns: {', '.join(missing_columns)}")

    df[args.ani_column] = pd.to_numeric(df[args.ani_column], errors="coerce")
    df[args.cov_column] = pd.to_numeric(df[args.cov_column], errors="coerce")
    df["__genome_id"] = (
        df["Genome_file"]
        .fillna("")
        .str.extract(r"/([^/]+)_genomic\.fna\.gz$", expand=False)
        .fillna(df["Genome_file"].fillna("").map(lambda value: Path(value).name if value else ""))
    )

    pass_hits = df.loc[(df[args.ani_column] >= args.ani) & (df[args.cov_column] >= args.cov)].copy()

    args.out_report.parent.mkdir(parents=True, exist_ok=True)
    pass_hits.drop(columns="__genome_id").to_csv(args.out_report, sep="\t", index=False)

    if args.out_references:
        references = sorted(pass_hits["Genome_file"].dropna().unique())
        args.out_references.parent.mkdir(parents=True, exist_ok=True)
        args.out_references.write_text("".join(f"{reference}\n" for reference in references))

    if args.out_summary:
        species_by_genome_id = load_species_by_genome_id(args.taxonomy_data)
        removed_genome_ids = sorted(set(df["__genome_id"].dropna()) - set(pass_hits["__genome_id"].dropna()))

        if removed_genome_ids:
            valid_removed_genome_ids = []
            for genome_id in removed_genome_ids:
                if genome_id:
                    valid_removed_genome_ids.append(genome_id)

            # Start a DataFrame with one row per removed genome ID.
            removed_genomes_df = pd.DataFrame({"genome_id": valid_removed_genome_ids})

            # Add the species name for each genome ID
            removed_genomes_df["species_name"] = removed_genomes_df["genome_id"].map(species_by_genome_id)
            removed_genomes_df["species_name"] = removed_genomes_df["species_name"].fillna("unknown_species")

            # Collapse removed genome IDs into a comma-separated list for each species.
            summary_df = (
                removed_genomes_df.groupby("species_name", sort=True)["genome_id"]
                .agg(lambda genome_ids: ",".join(sorted(set(genome_ids))))
                .reset_index(name="removed_genome_ids")
            )

            # Count how many unique removed references were recorded for each species.
            summary_df["removed_reference_count"] = summary_df["removed_genome_ids"].str.split(",").map(len)
            summary_df = summary_df[SUMMARY_COLUMNS]
        else:
            summary_df = pd.DataFrame(columns=SUMMARY_COLUMNS)

        args.out_summary.parent.mkdir(parents=True, exist_ok=True)
        summary_df.to_csv(args.out_summary, sep="\t", index=False)


if __name__ == "__main__":
    main()
