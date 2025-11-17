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
    except InvalidConfig:
        logging.error(
            f"Invalid config structure detected in configuration file {path}: {e}"
        )
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


def enrich_fields(df: pd.DataFrame) -> pd.DataFrame:
    df["sample_or_strain_name"] = (
        df[JOIN_COLUMN].str.rsplit("_", expand=True, n=1).loc[:, 0]
    )
    df["genome_status"] = df[JOIN_COLUMN].apply(
        lambda x: "mag" if "MAG" in x.upper() else "isolate"
    )
    return df


def merge_all(*dfs: pd.DataFrame) -> pd.DataFrame:
    df, other_dfs = dfs[0], dfs[1:]
    join_keys = [JOIN_COLUMN]
    for other_df in other_dfs:
        overlapping_cols = set(df.columns) & set(other_df.columns)
        cols_to_drop = overlapping_cols - set(join_keys)
        other_df_clean = other_df.drop(columns=cols_to_drop)
        df = df.merge(other_df_clean, on=join_keys, how="left")
    df = enrich_fields(df)
    return df


def reorder_columns(df: pd.DataFrame) -> pd.DataFrame:
    priority_cols = [
        JOIN_COLUMN,
        f"postqc_{JOIN_COLUMN_STEM}",
        "sample_or_strain_name",
        "genome_status",
    ]
    other_cols = sorted([col for col in df.columns if col not in priority_cols])
    new_order = priority_cols + other_cols
    return df[new_order]


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
        dfs.append(df)

    merged = merge_all(*dfs)
    merged = reorder_columns(merged)

    merged.to_csv(args.output, sep="\t", index=False)
    logging.info(f"Merged QC report written to: {args.output}")


if __name__ == "__main__":
    main()
