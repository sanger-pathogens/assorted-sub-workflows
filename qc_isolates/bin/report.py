#!/usr/bin/env python3

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Optional

import pandas as pd

JOIN_COLUMN = "genome_name"

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
    valid_tools = {"GUNC", "CHECKM2", "GTDBTK", "QUAST_SUMMARY"}
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


def process_input_report(df: pd.DataFrame, config: dict, tool: str) -> pd.DataFrame:
    tool = tool.lower()
    columns = config[tool.upper()]["keep_columns"]
    rename_dict = {}

    rename_dict[config[tool.upper()]["id_column"]] = JOIN_COLUMN

    for column in columns:
        rename_dict[column] = f"{tool}_{column}"

    return df.rename(columns=rename_dict)[rename_dict.values()]



def enrich_fields(df: pd.DataFrame) -> pd.DataFrame:
    df["sample_or_strain_name"] = (
        df[JOIN_COLUMN].str.rsplit("_", expand=True, n=1).loc[:, 0]
    )
    df["genome_status"] = "isolate"
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
        "sample_or_strain_name",
        "genome_status",
    ]
    other_cols = sorted([col for col in df.columns if col not in priority_cols])
    new_order = priority_cols + other_cols
    return df[new_order]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge reports from different QC tools and metadata."
    )
    # Args for tools must be consistently lowercase versions of the tool name in the config for main logic
    parser.add_argument(
        "--checkm2",
        type=Path,
        required=True,
        help="Path to CheckM2 TSV",
    )
    parser.add_argument(
        "--gunc", 
        type=Path, 
        required=True, 
        help="Path to GUNC TSV"
    )
    parser.add_argument(
        "--gtdbtk", 
        type=Path, 
        required=True, 
        help="Path to GTDBTk TSV"
    )
    parser.add_argument(
        "--quast", 
        type=Path, 
        required=True, 
        help="Path to QUAST_SUMMARY TSV"
    )
    parser.add_argument(
        "--config", 
        type=Path, 
        required=True, 
        help="Path to report configuration JSON file"
    )
    parser.add_argument(
        "--output", 
        type=Path, 
        required=True, 
        help="Output TSV path"
    )
    return parser.parse_args()


def main():
    setup_logging()
    args = parse_args()
    config = read_config(args.config)


    dfs = []
    
    # Processing in the order supplied in report config json
    # Order matters where two reports have a column of the same name, first is kept
    for tool in config.keys():
        tool = tool.lower()
        tsv_path = vars(args).get(tool)
        if tsv_path:
            df = process_input_report(read_tsv(tsv_path, tool), config, tool)
            dfs.append(df)

    merged = merge_all(*dfs)
    merged = reorder_columns(merged)

    merged.to_csv(args.output, sep="\t", index=False)
    logging.info(f"Merged QC report written to: {args.output}")


if __name__ == "__main__":
    main()
