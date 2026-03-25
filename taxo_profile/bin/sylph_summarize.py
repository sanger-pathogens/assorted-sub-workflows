#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd


def load_report(path):
    df = pd.read_csv(path, sep="\t")
    df["__report_path"] = str(path)
    # Sample ID is derived from the report filename prefix.
    df["__sample_id"] = Path(path).name.replace("_sylph_profile.tsv", "")
    return df


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--reports", nargs="+", required=True)
    p.add_argument("--ani", type=float, required=True)
    p.add_argument("--cov", type=float, required=True)
    p.add_argument("--ani-column", default="Adjusted_ANI")
    p.add_argument("--cov-column", default="Eff_cov")
    p.add_argument("--genome_path_prefix", default=None)
    p.add_argument("--out-references", default="references.txt")
    p.add_argument("--out-report", default="sylph_report.txt")
    p.add_argument("--out-summary", default="sylph_summary.tsv")
    args = p.parse_args()

    frames = [load_report(p) for p in args.reports]
    if not frames:
        raise SystemExit("No Sylph reports provided.")

    df = pd.concat(frames, ignore_index=True)

    ani_col = args.ani_column
    cov_col = args.cov_column

    # TODO: consider abundance thresholds once defined (e.g. Taxonomic_abundance or Sequence_abundance)

    pass_mask = (df[ani_col] >= args.ani) & (df[cov_col] >= args.cov)
    pass_hits = df.loc[pass_mask].copy()

    if args.genome_path_prefix:
        pass_hits["Genome_file"] = f"{Path(args.genome_path_prefix)}/" + pass_hits["Genome_file"]

    if len(pass_hits) > 0:
        pass_hits.to_csv(args.out_report, sep="\t", index=False)

    references = sorted(pass_hits["Genome_file"].dropna().unique())

    if references:
        Path(args.out_references).write_text("\n".join(references) + "\n")

    # Summarize per reference genome across all samples.
    summaries = []
    for genome_file, gdf in df.groupby("Genome_file", sort=True):
        gdf_pass = gdf[(gdf[ani_col] >= args.ani) & (gdf[cov_col] >= args.cov)]
        passed = int(not gdf_pass.empty)
        failed_thresholds = int(passed == 0)

        summaries.append(
            {
                "reference_genome": genome_file,
                "pass": passed,
                "failed_thresholds": failed_thresholds,
            }
        )

    summary_df = pd.DataFrame(summaries)

    # Append a TOTAL row to make quick tallying easy for reviews.
    totals = {
        "reference_genome": "TOTAL",
        "pass": int(summary_df["pass"].sum()),
        "failed_thresholds": int(summary_df["failed_thresholds"].sum()),
    }

    summary_df = pd.concat([summary_df, pd.DataFrame([totals])], ignore_index=True)
    summary_df.to_csv(args.out_summary, sep="\t", index=False)


if __name__ == "__main__":
    main()
