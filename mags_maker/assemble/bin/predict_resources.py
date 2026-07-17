#!/usr/bin/env python3

import argparse
import gzip
from pathlib import Path

import joblib
import pandas as pd


def count_reads(fastq):
    """Return number of reads in a FASTQ(.gz)."""
    opener = gzip.open if str(fastq).endswith(".gz") else open

    with opener(fastq, "rt") as fh:
        lines = sum(1 for _ in fh)

    return lines // 4


def load_richness(report):
    """
    Read richness from a Sylph report.
    """
    df = pd.read_csv(report, sep="\t")

    if "Richness" in df.columns:
        return float(df["Richness"].iloc[0])

    raise ValueError("Could not find Richness column in Sylph report.")


def memory_to_label(memory_mb):
    """
    Convert a predicted memory requirement (MB) to the smallest
    Nextflow memory label that satisfies it.
    """

    label_dict = {
        50: "mem_50M",
        100: "mem_100M",
        250: "mem_250M",
        500: "mem_500M",
        1024: "mem_1",
        2048: "mem_2",
        4096: "mem_4",
        8192: "mem_8",
        10240: "mem_10",
        16384: "mem_16",
        20480: "mem_20",
        32768: "mem_32",
        65536: "mem_64",
        98304: "mem_96",
        122880: "mem_120",
    }

    for threshold, label in label_dict.items():
        if memory_mb <= threshold:
            return label

    # Anything larger than largest configured label
    return "mem_120"


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--model", required=True)
    parser.add_argument("--fastq", required=True)
    parser.add_argument("--sylph", required=True)

    args = parser.parse_args()

    model = joblib.load(args.model)

    reads = count_reads(args.fastq)
    richness = load_richness(args.sylph)

    X = pd.DataFrame(
        {
            "total_reads": [reads],
            "Richness": [richness],
        }
    )

    predicted_memory = float(model.predict(X)[0])

    label = memory_to_label(predicted_memory)

    print(label)


if __name__ == "__main__":
    main()