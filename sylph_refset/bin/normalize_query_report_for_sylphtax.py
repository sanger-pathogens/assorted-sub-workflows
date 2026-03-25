#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True, help="Path to raw sylph query report TSV.")
    parser.add_argument("--output", type=Path, required=True, help="Path to normalized TSV for sylph-tax taxprof.")
    return parser.parse_args()


def main():
    args = parse_args()
    columns = [
        "Sample_file",
        "Genome_file",
        "Taxonomic_abundance",
        "Sequence_abundance",
        "Adjusted_ANI",
        "Eff_cov",
        "ANI_5-95_percentile",
        "Eff_lambda",
        "Lambda_5-95_percentile",
        "Median_cov",
        "Mean_cov_geq1",
        "Containment_ind",
        "Naive_ANI",
        "kmers_reassigned",
        "Contig_name",
    ]

    with args.input.open(newline="") as src, args.output.open("w", newline="") as dst:
        reader = csv.DictReader(src, delimiter="\t")
        writer = csv.DictWriter(dst, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in reader:
            writer.writerow({col: row.get(col, "") for col in columns})


if __name__ == "__main__":
    main()
