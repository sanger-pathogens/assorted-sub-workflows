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
from pathlib import Path
import sys
import logging

import pod5 as p5

logger = logging.getLogger(__name__)

DEFAULT_VERSION = "1.3.1"
SCHEMA_VERSION = "1.0"

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

def classify(run_info):
    for rule in CHEMISTRY_RULES:
        if rule["match"](run_info):
            return rule["name"], rule["dorado_version"], rule["reason"]
    return "standard", DEFAULT_VERSION, ""

def build_report(pod5_file: Path, read) -> dict:
    run_info = read.run_info
    rule_name, dorado_version, reason = classify(run_info)
    return {
        "schema_version": SCHEMA_VERSION,
        "pod5_file": str(pod5_file),
        "rule": rule_name,
        "dorado_version": dorado_version,
        "sample_rate": getattr(run_info, "sample_rate", None),
        "flow_cell_product_code": getattr(run_info, "flow_cell_product_code", None),
        "sequencing_kit": getattr(run_info, "sequencing_kit", None),
        "reason": reason,
    }

def main() -> int:
    parser = argparse.ArgumentParser(description="Inspect the first read of a pod5 file.")
    parser.add_argument("pod5_file", type=Path, help="Path to a single pod5 file to inspect")
    parser.add_argument(
        "-o", "--output",
        type=Path,
        default=Path("read_info.json"),
        help="Path to write the JSON report to (default: %(default)s)",
    )
    parser.add_argument("-q", "--quiet", action="store_true", help="Suppress stdout report")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    if not args.pod5_file.is_file():
        logger.error("%s does not exist or is not a file", args.pod5_file)
        return 1

    try:
        with p5.Reader(args.pod5_file) as reader:
            read = next(reader.reads(), None)
            if read is None:
                logger.error("%s contains no reads", args.pod5_file)
                return 1
            report = build_report(args.pod5_file, read)
    except Exception:
        logger.exception("Failed to read %s", args.pod5_file)
        return 1


    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as fh:
        json.dump(report, fh, indent=2)

    if not args.quiet:
        print(json.dumps(report, indent=2))

    return 0


if __name__ == "__main__":
    sys.exit(main())