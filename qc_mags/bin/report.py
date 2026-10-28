#!/usr/bin/env python3
"""Merge the pre-/post-QC QUAST, CheckM2 and GUNC reports and the GTDB-Tk summary of one batch into one table."""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Optional

import pandas as pd

JOIN_COLUMN_STEM = "genome_name"
JOIN_COLUMN = f"preqc_{JOIN_COLUMN_STEM}"
POSTQC_COLUMN = f"postqc_{JOIN_COLUMN_STEM}"
TOOLS = ("GUNC", "CHECKM2", "GTDBTK", "QUAST")
COMPLETENESS_COLUMN = "checkm2_postqc_Completeness"
CONTAMINATION_COLUMN = "checkm2_postqc_Contamination"
GUNC_COLUMN = "gunc_postqc_pass.GUNC"
QUALITY_COLUMN = "mimag_quality"
REVIEW_COLUMN = "qc_review"
TAXONOMY_COLUMN = "gtdbtk_classification"
ISOLATE_REVIEW_FLAG = "isolate_completeness_le_90: check CheckM2 model coverage (novel species?) before trusting the tier"


def read_config(path: Path) -> dict:
    with open(path) as f:
        config = json.load(f)
    missing_tools = set(TOOLS) - config.keys()
    if missing_tools:
        raise SystemExit(f"Report config {path} is missing tool(s) {sorted(missing_tools)}")
    for tool in TOOLS:
        missing_keys = {"id_column", "keep_columns"} - config[tool].keys()
        if missing_keys:
            raise SystemExit(f"Report config {path}: tool '{tool}' is missing {sorted(missing_keys)}")
    unexpected = config.keys() - set(TOOLS)
    if unexpected:
        logging.warning(f"Report config contains keys that were not expected (and are ignored): {unexpected}")
    return config


def strip_postqc_name(names: pd.Series) -> pd.Series:
    """cleaned_<name>_filtered_kept_contigs -> <name>; names without the wrapper are left unchanged."""
    return names.str.extract(r"^cleaned_(.*)_filtered_kept_contigs$")[0].fillna(names)


def process_input_report(df: pd.DataFrame, tool_config: dict, prefix: str, id_column: str) -> pd.DataFrame:
    rename = {tool_config["id_column"]: id_column}
    rename.update({column: f"{prefix}_{column}" for column in tool_config["keep_columns"]})
    df = df.rename(columns=rename)[list(rename.values())]
    if id_column == POSTQC_COLUMN:
        df[JOIN_COLUMN] = strip_postqc_name(df[id_column])
    return df


def mimag_quality(completeness, contamination, gunc_pass, genome_type) -> Optional[str]:
    """Quality tier from post-QC CheckM2 plus post-QC GUNC (the rRNA/tRNA criteria
    of the full MIMAG HQ standard are not assessed). Same CheckM2 thresholds for
    every genome type:
      HQ: completeness > 90 and contamination < 5;
      MQ: completeness >= 50 and contamination < 10; otherwise LQ.
    GUNC (binning chimerism) gates MAGs only: HQ needs an explicit GUNC pass, a
    GUNC fail is LQ, no GUNC result is MQ at most. For isolates (pure cultures:
    no binning, so no binning chimeras) GUNC is advisory and never changes the tier.
    Empty if CheckM2 is missing."""
    if pd.isna(completeness) or pd.isna(contamination):
        return None
    isolate = genome_type == "isolate"
    g = "" if gunc_pass is None or pd.isna(gunc_pass) else str(gunc_pass).strip().lower()
    if not isolate and g == "false":
        return "LQ"
    if completeness > 90 and contamination < 5 and (isolate or g == "true"):
        return "HQ"
    if completeness >= 50 and contamination < 10:
        return "MQ"
    return "LQ"


def qc_review_flag(completeness, genome_type) -> str:
    """Isolates at completeness <= 90 are MQ, not HQ. CheckM2 is weakest on rare and novel species, so such an
    isolate may be a complete genome that the models under-score: flag it for review. The tier is not changed."""
    if genome_type != "isolate" or pd.isna(completeness):
        return ""
    return ISOLATE_REVIEW_FLAG if completeness <= 90 else ""


def add_mimag_quality(df: pd.DataFrame) -> pd.DataFrame:
    if COMPLETENESS_COLUMN not in df.columns or CONTAMINATION_COLUMN not in df.columns:
        logging.warning(f"{COMPLETENESS_COLUMN}/{CONTAMINATION_COLUMN} not in report; {QUALITY_COLUMN} left empty")
        df[QUALITY_COLUMN] = None
        df[REVIEW_COLUMN] = ""
        return df
    comp = pd.to_numeric(df[COMPLETENESS_COLUMN], errors="coerce")
    cont = pd.to_numeric(df[CONTAMINATION_COLUMN], errors="coerce")
    if GUNC_COLUMN in df.columns:
        gunc = df[GUNC_COLUMN]
    else:
        logging.warning(f"{GUNC_COLUMN} not in report; no MAG can be HQ (needs a GUNC pass)")
        gunc = [None] * len(df)
    gtype = df["genome_status"]
    df[QUALITY_COLUMN] = [mimag_quality(a, b, g, t) for a, b, g, t in zip(comp, cont, gunc, gtype)]
    df[REVIEW_COLUMN] = [qc_review_flag(a, t) for a, t in zip(comp, gtype)]
    return df


def read_genome_types(path) -> dict:
    """Per-genome types for mixed runs (columns genome, genome_type); validated by the workflow (genomeTypeMap)."""
    if path is None:
        return {}
    t = pd.read_csv(path, sep="\t", dtype=str)
    return dict(zip(t["genome"], t["genome_type"].str.strip().str.lower()))


def enrich_fields(df: pd.DataFrame, genome_type: str, genome_types: dict) -> pd.DataFrame:
    df["sample_or_strain_name"] = df[JOIN_COLUMN].str.rsplit("_", expand=True, n=1).loc[:, 0]
    df["genome_status"] = genome_type
    if genome_types:
        # mixed run: per-genome type from the table; genomes not listed keep the run-level type
        hit = df[JOIN_COLUMN].map(genome_types)
        df.loc[hit.notna(), "genome_status"] = hit[hit.notna()]
        logging.info(f"--genome_types: {int(hit.notna().sum())} of {len(df)} genomes typed from the table "
                     f"({int((hit == 'isolate').sum())} isolate); the rest use --genome_type {genome_type}")
        if hit.notna().sum() == 0:
            logging.warning("--genome_types matched no genome in this report: check the 'genome' column "
                            "uses input genome names without extension (preqc_genome_name)")
    return add_mimag_quality(df)


def merge_all(dfs: list) -> pd.DataFrame:
    df = dfs[0]
    for other in dfs[1:]:
        # post-QC reports all carry postqc_genome_name; keep only the first copy
        df = df.merge(other.drop(columns=(set(df.columns) & set(other.columns)) - {JOIN_COLUMN}),
                      on=JOIN_COLUMN, how="left")
    return df


def reorder_columns(df: pd.DataFrame) -> pd.DataFrame:
    priority_cols = [JOIN_COLUMN, POSTQC_COLUMN, "sample_or_strain_name", "genome_status", QUALITY_COLUMN,
                     REVIEW_COLUMN]
    return df[priority_cols + sorted(col for col in df.columns if col not in priority_cols)]


def check_taxonomy(df: pd.DataFrame) -> None:
    """Warn about genomes without a GTDB-Tk call; fail if NO genome has one, as that means the GTDB-Tk summary
    is empty or did not join on genome name, and the report would look complete without any taxonomy."""
    if df.empty:
        logging.warning("Merged report has no genomes; skipping taxonomy check")
        return
    if TAXONOMY_COLUMN in df.columns:
        missing = df[TAXONOMY_COLUMN].isna() | (df[TAXONOMY_COLUMN].astype(str).str.strip() == "")
    else:
        missing = pd.Series(True, index=df.index)
    n_missing = int(missing.sum())
    if n_missing:
        logging.warning(
            f"{n_missing}/{len(df)} genomes lack {TAXONOMY_COLUMN}: "
            + ", ".join(df.loc[missing, JOIN_COLUMN].astype(str).head(20))
            + (" ..." if n_missing > 20 else "")
        )
    if n_missing == len(df):
        logging.error(f"ALL {len(df)} genomes lack {TAXONOMY_COLUMN}: the GTDB-Tk summary is empty or did not "
                      "join on genome name. Refusing to write a report without taxonomy.")
        sys.exit(2)


# argument -> (tool, column prefix, name of the id column after renaming), in merge order (pre-QC QUAST first)
REPORTS = {
    "pre_qc_quast": ("QUAST", "quast_preqc", JOIN_COLUMN),
    "pre_qc_checkm2": ("CHECKM2", "checkm2_preqc", JOIN_COLUMN),
    "pre_qc_gunc": ("GUNC", "gunc_preqc", JOIN_COLUMN),
    "gtdbtk": ("GTDBTK", "gtdbtk", JOIN_COLUMN),
    "post_qc_quast": ("QUAST", "quast_postqc", POSTQC_COLUMN),
    "post_qc_checkm2": ("CHECKM2", "checkm2_postqc", POSTQC_COLUMN),
    "post_qc_gunc": ("GUNC", "gunc_postqc", POSTQC_COLUMN),
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in REPORTS:
        parser.add_argument(f"--{arg}", type=Path, required=True, help=f"{arg} report TSV")
    parser.add_argument("--config", type=Path, required=True, help="Path to report configuration file")
    parser.add_argument("--output", type=Path, required=True, help="Output TSV path")
    parser.add_argument("--genome_type", choices=("mag", "isolate"), default="mag",
                        help="genome type of every genome in the run: 'mag' (default; GUNC gates the tier) or "
                        "'isolate' (pure cultures; GUNC reported but advisory)")
    parser.add_argument("--genome_types", type=Path, default=None,
                        help="optional TSV for mixed runs: columns genome (input genome name without extension) "
                        "and genome_type (mag|isolate); genomes not listed use --genome_type")
    return parser.parse_args()


def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    args = parse_args()
    config = read_config(args.config)

    dfs = []
    for arg, (tool, prefix, id_column) in REPORTS.items():
        df = process_input_report(pd.read_csv(getattr(args, arg), sep="\t"), config[tool], prefix, id_column)
        if tool == "GTDBTK":
            # GTDB-Tk ids are post-QC names when it runs on the cleaned fastas; strip them to the pre-QC join key
            df[JOIN_COLUMN] = strip_postqc_name(df[JOIN_COLUMN])
        dfs.append(df)

    merged = merge_all(dfs)
    merged = enrich_fields(merged, args.genome_type, read_genome_types(args.genome_types))
    merged = reorder_columns(merged)
    check_taxonomy(merged)

    merged.to_csv(args.output, sep="\t", index=False)
    logging.info(f"Merged QC report written to: {args.output}")


if __name__ == "__main__":
    main()
