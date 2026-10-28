#!/usr/bin/env python3

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Optional

import pandas as pd

JOIN_COLUMN_STEM = "genome_name"
JOIN_COLUMN = f"preqc_{JOIN_COLUMN_STEM}"
COMPLETENESS_COLUMN = "checkm2_postqc_Completeness"
CONTAMINATION_COLUMN = "checkm2_postqc_Contamination"
GUNC_COLUMN = "gunc_postqc_pass.GUNC"
QUALITY_COLUMN = "mimag_quality"
REVIEW_COLUMN = "qc_review"
TAXONOMY_COLUMN = "gtdbtk_classification"
ISOLATE_REVIEW_FLAG = "isolate_completeness_le_90: check CheckM2 model coverage (novel species?) before trusting the tier"


def setup_logging(log_file: str = "qc_merge.log"):
    logging.basicConfig(
        level=logging.INFO,
        handlers=[logging.StreamHandler(), logging.FileHandler(log_file, mode="w")],
        format="%(asctime)s - %(levelname)s - %(message)s",
    )
    logging.info("Logging initialized.")


def read_tsv(path: Path, name: str) -> pd.DataFrame:
    try:
        df = pd.read_csv(path, sep="\t")
        logging.info(f"{name} loaded successfully: {path}")
    except Exception as e:
        logging.error(f"Error reading {name} from {path}: {e}")
        sys.exit(1)
    return df


def read_config(path: Path) -> dict:
    try:
        with open(path) as f:
            config = json.load(f)
    except Exception as e:
        logging.error(f"Error reading configuration file {path}: {e}")
        sys.exit(1)
    try:
        validate_config(config)
    except InvalidConfig as e:
        logging.error(
            f"Invalid config structure detected in configuration file {path}: {e}"
        )
        sys.exit(1)
    return config


class InvalidConfig(ValueError):
    pass


def validate_config(config: dict) -> None:
    valid_tools = {"GUNC", "CHECKM2", "GTDBTK", "QUAST"}
    expected_tools = valid_tools & config.keys()
    unexpected_tools = set(config.keys()) - valid_tools
    if expected_tools == set():
        raise InvalidConfig(
            f"Config does not contain at least one valid tool: {valid_tools}"
        )
    if unexpected_tools != set():
        logging.warning(
            f"Config contains keys that were not expected (and are ignored): {unexpected_tools}"
        )
    for tool in expected_tools:
        missing_keys = {"id_column", "keep_columns"} - config[tool].keys()
        if missing_keys != set():
            raise InvalidConfig(
                f"Config is missing a mandatory key for tool '{tool}': {missing_keys}"
            )
    logging.info("Config passed validation")


def process_input_report(
    df: pd.DataFrame, config: dict, tool: str, qc_stage: Optional[str] = None
) -> pd.DataFrame:
    tool = tool.lower()
    columns = config[tool.upper()]["keep_columns"]
    rename_dict = {}
    if qc_stage is not None:
        assert qc_stage in ("preqc", "postqc")
        rename_dict[config[tool.upper()]["id_column"]] = (
            f"{qc_stage}_{JOIN_COLUMN_STEM}"
        )
    else:
        rename_dict[config[tool.upper()]["id_column"]] = JOIN_COLUMN
    for column in columns:
        if qc_stage is not None:
            rename_dict[column] = f"{tool}_{qc_stage}_{column}"
        else:
            rename_dict[column] = f"{tool}_{column}"
    return df.rename(columns=rename_dict)[rename_dict.values()]


def process_postqc_genome_name(df: pd.DataFrame) -> pd.DataFrame:
    new_df = df.copy()
    new_df[JOIN_COLUMN] = df[f"postqc_{JOIN_COLUMN_STEM}"].str.extract(
        r"cleaned_(.*)_filtered_kept_contigs"
    )
    return new_df


def enrich_fields(df: pd.DataFrame, genome_type: str, genome_types: dict) -> pd.DataFrame:
    df["sample_or_strain_name"] = (
        df[JOIN_COLUMN].str.rsplit("_", expand=True, n=1).loc[:, 0]
    )
    # genome type comes from --genome_type/--genome_types (curated provenance), never from the genome name
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


def merge_all(*dfs: pd.DataFrame, genome_type: str = "mag", genome_types: dict = None) -> pd.DataFrame:
    df, other_dfs = dfs[0], dfs[1:]
    join_keys = [JOIN_COLUMN]
    for other_df in other_dfs:
        overlapping_cols = set(df.columns) & set(other_df.columns)
        cols_to_drop = overlapping_cols - set(join_keys)
        other_df_clean = other_df.drop(columns=cols_to_drop)
        df = df.merge(other_df_clean, on=join_keys, how="left")
    df = enrich_fields(df, genome_type, genome_types)
    return df


def reorder_columns(df: pd.DataFrame) -> pd.DataFrame:
    priority_cols = [
        JOIN_COLUMN,
        f"postqc_{JOIN_COLUMN_STEM}",
        "sample_or_strain_name",
        "genome_status",
        QUALITY_COLUMN,
        REVIEW_COLUMN,
    ]
    other_cols = sorted([col for col in df.columns if col not in priority_cols])
    new_order = priority_cols + other_cols
    return df[new_order]


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


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge GUNC, CheckM2 and post-QC metadata."
    )
    parser.add_argument(
        "--pre_qc_quast", type=Path, required=True, help="Path to pre-QC QUAST report TSV"
    )
    parser.add_argument(
        "--pre_qc_checkm2", type=Path, required=True, help="Path to pre-QC CheckM2 TSV"
    )
    parser.add_argument(
        "--pre_qc_gunc", type=Path, required=True, help="Path to pre-QC GUNC TSV"
    )
    parser.add_argument(
        "--post_qc_quast", type=Path, required=True, help="Path to post-QC QUAST report TSV"
    )
    parser.add_argument(
        "--post_qc_checkm2", type=Path, required=True, help="Path to post-QC CheckM2 TSV",
    )
    parser.add_argument(
        "--post_qc_gunc", type=Path, required=True, help="Path to post-QC GUNC TSV"
    )
    parser.add_argument(
        "--gtdbtk", type=Path, required=True, help="Path to GTDBTK TSV"
    )
    parser.add_argument(
        "--config", type=Path, required=True, help="Path to report configuration file"
    )
    parser.add_argument(
        "--output", type=Path, required=True, help="Output TSV path"
    )
    parser.add_argument(
        "--genome_type", choices=("mag", "isolate"), default="mag",
        help="genome type of every genome in the run: 'mag' (default; GUNC gates the tier) or "
        "'isolate' (pure cultures; GUNC reported but advisory)",
    )
    parser.add_argument(
        "--genome_types", type=Path, default=None,
        help="optional TSV for mixed runs: columns genome (input genome name without extension) and "
        "genome_type (mag|isolate); genomes not listed use --genome_type",
    )
    return parser.parse_args()


def main():
    setup_logging()
    args = parse_args()

    config = read_config(args.config)

    args_qc_stage = {
        "pre_qc_quast": ["preqc","quast"],
        "pre_qc_checkm2": ["preqc", "checkm2"],
        "pre_qc_gunc": ["preqc", "gunc"],
        "gtdbtk": [None, "gtdbtk"],
        "post_qc_quast": ["postqc", "quast"],
        "post_qc_checkm2": ["postqc", "checkm2"],
        "post_qc_gunc": ["postqc", "gunc"]
    }

    dfs = []
    for arg, qc in args_qc_stage.items():
        qc_stage = qc[0]
        tool = qc[1]
        df = read_tsv(vars(args)[arg], arg)
        df = process_input_report(df, config, tool, qc_stage)

        if qc_stage == "postqc":
            df = process_postqc_genome_name(df)
        elif tool == "gtdbtk":
            # GTDB-Tk ids are post-QC names when it runs on the cleaned fastas; strip them to the pre-QC join key
            stripped = df[JOIN_COLUMN].str.extract(r"^cleaned_(.*)_filtered_kept_contigs$").loc[:, 0]
            df[JOIN_COLUMN] = stripped.fillna(df[JOIN_COLUMN])
        dfs.append(df)

    merged = merge_all(*dfs, genome_type=args.genome_type, genome_types=read_genome_types(args.genome_types))
    merged = reorder_columns(merged)

    check_taxonomy(merged)

    merged.to_csv(args.output, sep="\t", index=False)
    logging.info(f"Merged QC report written to: {args.output}")


if __name__ == "__main__":
    main()
