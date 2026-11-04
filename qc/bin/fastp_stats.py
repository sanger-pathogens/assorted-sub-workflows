#!/usr/bin/env python3
"""Collate fastp JSON reports (<ID>.fastp.json) into one TSV, one row per sample, written to stdout."""
import argparse
import json
from pathlib import Path

COLUMNS = [
    "Sample_id",
    "reads_before",
    "reads_after",
    "pct_reads_kept",
    "bases_before",
    "bases_after",
    "q30_rate_after",
    "gc_content_after",
    "duplication_rate",
    "low_quality_reads",
    "too_many_N_reads",
    "too_short_reads",
    "adapter_trimmed_reads",
]


def sample_row(json_path):
    report = json.loads(Path(json_path).read_text())
    before = report["summary"]["before_filtering"]
    after = report["summary"]["after_filtering"]
    filtering = report["filtering_result"]
    kept = 100 * after["total_reads"] / before["total_reads"] if before["total_reads"] else 0
    return [
        Path(json_path).name.removesuffix(".fastp.json"),
        before["total_reads"],
        after["total_reads"],
        f"{kept:.2f}",
        before["total_bases"],
        after["total_bases"],
        after["q30_rate"],
        after["gc_content"],
        report.get("duplication", {}).get("rate", ""),
        filtering["low_quality_reads"],
        filtering["too_many_N_reads"],
        filtering["too_short_reads"],
        report.get("adapter_cutting", {}).get("adapter_trimmed_reads", 0),
    ]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fastp_jsons", nargs="+", help="fastp JSON reports named <ID>.fastp.json")
    args = parser.parse_args()

    print("\t".join(COLUMNS))
    for row in sorted(sample_row(path) for path in args.fastp_jsons):
        print("\t".join(str(value) for value in row))


if __name__ == "__main__":
    main()
