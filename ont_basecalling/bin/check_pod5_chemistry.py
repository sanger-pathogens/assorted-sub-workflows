#!/usr/bin/env python3
"""
check_pod5_chemistry.py

Take the first read and check its run_info against a curated list of chemistry rules to determine
whether the pod5 file is compatible with the current Dorado release. reports as json for nextflow

This is a curated exception list of chemistries we've *observed* to fail
with current Dorado releases

Requires manual Update CHEMISTRY_RULES reactively as new failures are found.
"""
import argparse
import json
import sys

import pod5 as p5

#Rules for pod5s chemistry compatibility with dorado. If a pod5 matches a rule, the rule name and dorado version 
#will be reported in the output json otherwise the pod5 is assumed to be compatible with the current dorado release
CHEMISTRY_RULES = [
    {
        "name": "r10.4.1_4khz",
        "match": lambda ri: (
            ri.sample_rate == 4000
            and (ri.flow_cell_product_code or "").upper().startswith("FLO-MIN114")
        ),
        "dorado_version": "0.9.6",
        "reason": "DNA r10.4.1 e8.2 4kHz chemistry deprecated since Dorado v1.0.0 (the 4000Hz sample rate is no longer supported)",
    },
    {
        "name": "r9.4.1_any_kit",
        "match": lambda ri: (ri.flow_cell_product_code or "").upper().startswith("FLO-MIN106"),
        "dorado_version": "0.9.6",
        "reason": "R9.4.1 models removed from newer Dorado releases",
    },
]

DEFAULT_VERSION = "current"


def classify(run_info):
    for rule in CHEMISTRY_RULES:
        if rule["match"](run_info):
            return rule["name"], rule["dorado_version"], rule["reason"]
    return "standard", DEFAULT_VERSION, ""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("pod5_file", help="Path to a single pod5 file to inspect")
    parser.add_argument(
        "-o", "--output",
        default="read_info.json",
        help="Path to write the JSON report to (default: %(default)s)",
    )
    args = parser.parse_args()

    with p5.Reader(args.pod5_file) as reader:
        try:
            read = next(reader.reads())
        except StopIteration:
            sys.exit(f"ERROR: {args.pod5_file} contains no reads")

        run_info = read.run_info
        rule_name, dorado_version, reason = classify(run_info)

        report = {
            "pod5_file": str(args.pod5_file),
            "rule": rule_name,
            "dorado_version": dorado_version,
            "sample_rate": run_info.sample_rate,
            "flow_cell_product_code": run_info.flow_cell_product_code,
            "sequencing_kit": run_info.sequencing_kit,
            "reason": reason,
        }

    with open(args.output, "w") as fh:
        json.dump(report, fh, indent=2)

    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()