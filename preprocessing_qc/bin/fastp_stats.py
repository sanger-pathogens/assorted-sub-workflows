#!/usr/bin/env python3
"""Collate fastp JSON reports (<ID>.fastp.json) into one TSV, one row per sample, written to stdout."""
import argparse
import json
from pathlib import Path

COLUMNS_STATS = [
    "sample_id",
    "reads_before",
    "reads_after",
    "pct_reads_kept",
    "bases_before",
    "bases_after",
    "q30_rate_before",
    "q30_rate_after",
    "gc_content_before",
    "gc_content_after",
    "duplication_rate",
    "low_quality_reads_filtered",
    "too_many_N_reads_filtered",
    "too_short_reads_filtered",
    "adapter_trimmed_reads",
]

COLUMNS_FAIL_PASS = [
  "sample_id",
  "read1_per_base_quality",
  "read1_per_seq_quality",
  "read1_per_base_N",
  "read2_per_base_quality",
  "read2_per_seq_quality",
  "read2_per_base_N",
]

def fail_pass_fastp_filtered_reads(report, sample_id, outfailpasshandle):
    grades = [sample_id]
    for read in ("read1_after_filtering", "read2_after_filtering"):
        if read not in report:
            grades.extend(["NA", "NA", "NA"])
            continue
        r = report[read]
        qual = r["quality_curves"]["mean"]
        n = r["content_curves"]["N"]
        if not qual or not n:
            # no reads left after filtering
            grades.extend(["FAIL", "FAIL", "FAIL"])
            continue

        pos_qual = min(qual)                 # per base quality proxy
        mean_qual = sum(qual) / len(qual)    # per sequence quality proxy
        max_n = max(n)                       # per base N content

        def grade(v, warn, fail, higher_is_bad=False):
            if higher_is_bad:
                return "FAIL" if v > fail else "WARN" if v > warn else "PASS"
            return "FAIL" if v < fail else "WARN" if v < warn else "PASS"

        grades.append(grade(pos_qual, 25, 20))
        grades.append(grade(mean_qual, 27, 20))
        grades.append(grade(max_n, 0.05, 0.20, higher_is_bad=True))

    outfailpasshandle.write("\t".join(str(value) for value in grades) + "\n")

def extract_stats_to_json(report, json_path, outtsvhandle):
    before = report["summary"]["before_filtering"]
    after = report["summary"]["after_filtering"]
    filtering = report["filtering_result"]
    kept = 100 * after["total_reads"] / before["total_reads"] if before["total_reads"] else 0
    outtsvhandle.write("\t".join(str(value) for value in [
        Path(json_path).name.removesuffix(".fastp.json"),
        before["total_reads"],
        after["total_reads"],
        f"{kept:.2f}",
        before["total_bases"],
        after["total_bases"],
        before["q30_rate"],
        after["q30_rate"],
        before["gc_content"],
        after["gc_content"],
        report.get("duplication", {}).get("rate", ""),
        filtering["low_quality_reads"],
        filtering["too_many_N_reads"],
        filtering["too_short_reads"],
        report.get("adapter_cutting", {}).get("adapter_trimmed_reads", 0),
    ]) + "\n")

def sample_row(json_path, outtsvhandle, outfailpasshandle):
    report = json.loads(Path(json_path).read_text())
    sample_id = Path(json_path).name.removesuffix(".fastp.json")
    extract_stats_to_json(report, json_path, outtsvhandle)
    fail_pass_fastp_filtered_reads(report, sample_id, outfailpasshandle)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fastp_jsons", nargs="+", help="fastp JSON reports named <ID>.fastp.json")
    args = parser.parse_args()

    tsvout = open("fastp_statistics.tsv", "w")
    tsvout.write("\t".join(COLUMNS_STATS) + "\n")
    failpassout = open("fastp_fail_pass.tsv", "w")
    failpassout.write("\t".join(COLUMNS_FAIL_PASS) + "\n")
    for path in sorted(args.fastp_jsons, key=lambda p: Path(p).name):
        sample_row(path, tsvout, failpassout)
    tsvout.close()
    failpassout.close()


if __name__ == "__main__":
    main()
