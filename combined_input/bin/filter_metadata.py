#!/usr/bin/env python3

import argparse
import datetime as dt
import pandas as pd
import logging
import sys

from pandas.api.types import infer_dtype
from pathlib import Path


def parse_arguments():
    timestamp = dt.datetime.now().strftime("%Y-%M-%d_%H-%M-%S")

    parser = argparse.ArgumentParser(
        description="Filter rows of a TSV file based on conditions.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--input",
        "-i",
        type=Path,
        required=True,
        help="Path to the input TSV file."
    )
    parser.add_argument(
        "--filter_manifest",
        "-f",
        type=Path,
        required=True,
        help="Path to a TSV manifest specifying the columns in the input TSV to which filter and datatype conversion should be applied. The manifest should contain 3 columns: column, filter, datatype. Entries in the 'column' column should match a column in the input TSV and be unique. The filter should match a string that could be supplied to `pd.DataFrame.query()`, e.g. 'age > 30'. The datatype can be int, float, datetime, bool and str.",
    )
    parser.add_argument(
        "--select",
        "-s",
        type=str,
        nargs="+",
        help="Specify columns to select in the output DataFrame. By default, all columns will be selected.",
    )
    parser.add_argument(
        "--missing_values",
        "-m",
        type=str,
        nargs="+",
        help="Specify values that should be interpreted as missing values.",
    )
    parser.add_argument(
        "--remove_header",
        "-r",
        action="store_true",
        help="Remove the header from the output TSV.",
    )
    parser.add_argument(
        "--error_on_invalid_type",
        "-e",
        action="store_true",
        help="During type conversion, upon encountering a value in the column that cannot be converted to the given datatype, raise an error. Default behaviour is to remove the row that contains an invalid value.",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        help="Path to the output file to save the filtered TSV. If not specified, the metadata DataFrame will be printed to stdout.",
    )
    parser.add_argument(
        "--logfile",
        "-l",
        type=str,
        help="Path to the log file, defaults to filter_metadata-<timestamp>.log",
        default=f"filter_metadata-{timestamp}.log",
    )

    return parser.parse_args()


def build_query_string(conditions: list[str]) -> str:
    """Converts filter conditions into a valid query string for pandas"""
    query_list = []
    for condition in conditions:
        query_list.append(condition)
    return " and ".join(query_list)


def apply_filters(df: pd.DataFrame, filters: list[str]) -> pd.DataFrame:
    """Apply filters to select rows from the given DataFrame"""
    query_string = build_query_string(filters)
    try:
        filtered_df = df.query(query_string)
    except Exception as e:
        logging.error(f"Error while filtering with query '{query_string}': {e}")
        logging.info(f"Are you sure the columns are of correct dtype?")
        logging.info("Try forcing column dtype using the -d option.")
        sys.exit(1)
    return filtered_df


def parse_filter_manifest(manifest: Path, sep: str = "\t") -> pd.DataFrame:
    """
    Parses a manifest file specifying column, filter and data type into a pandas DataFrame.

    Parameters:
    manifest (Path): Path to the manifest file

    Returns:
    pd.DataFrame: A three column dataframe specifying column, filter and data type
    """
    if not manifest.is_file():
        logging.error(f"The given manfiest {manifest} is not a file.")
        sys.exit(1)

    df = pd.read_csv(manifest, sep=sep, skip_blank_lines=True)

    errors = 0

    # Check invalid headings
    headings = set(df.columns)
    missing_headings, unrecognized_headings = validate_headings(headings)
    if missing_headings:
        logging.error(
            f"Missing required headings from manifest {manifest}: {missing_headings}"
        )
        errors += 1
    if unrecognized_headings:
        logging.error(
            f"Unrecognized headings from manifest {manifest}: {unrecognized_headings}"
        )
        errors += 1

    # Check missing values
    num_missing_values = df.isnull().sum().sum()
    if num_missing_values != 0:
        logging.error(
            f"{num_missing_values} missing values were detected in the manifest {manifest}. Please ensure the manifest contains no missing values."
        )
        errors += 1

    # Check uniqueness of column filter specifications
    column_manifest_heading = "column"
    columns_names = df[column_manifest_heading]
    duplicate_column_names_mask = columns_names.duplicated(keep="first")
    if sum(duplicate_column_names_mask) != 0:
        logging.error(
            f"Found duplicate values in column {column_manifest_heading}:\n {df[duplicate_column_names_mask]}"
        )
        errors += 1

    if errors:
        logging.error(f"Found {errors} errors. Please check the manifest.")
        sys.exit(1)

    return df


def validate_headings(headings: set):
    valid_headings = {"column", "filter", "datatype"}
    missing_headings = valid_headings.difference(headings)
    unrecognized_headings = set(headings).difference(valid_headings)
    return missing_headings, unrecognized_headings


def safe_convert_column(
    df: pd.DataFrame, column: str, dtype: str, error_on_invalid_type: bool = False
) -> pd.DataFrame:
    """
    Attempts to convert a DataFrame column to the specified data type.
    Removes rows where any value prevents conversion.

    Parameters:
    df (pd.DataFrame): The DataFrame containing the column to convert.
    column (str): The name of the column to convert.
    dtype (str): The target data type (e.g., 'int', 'float', 'datetime', etc.).
    error_on_invalid (bool): If True, error on invalid values that cannot be converted, instead of row removal.

    Returns:
    pd.DataFrame: The DataFrame with the column converted, with rows with invalid values removed (if necessary).
    """
    errors_arg = "raise" if error_on_invalid_type else "coerce"
    try:
        if dtype == "int" or dtype == "float":
            # Use pd.to_numeric with errors='coerce' for numeric conversion
            df.loc[:, column] = pd.to_numeric(df[column], errors=errors_arg)
        elif dtype == "datetime":
            # Use pd.to_datetime with errors='coerce' for datetime conversion
            df[column] = pd.to_datetime(df[column], errors=errors_arg)
        else:
            # For other types, try using astype with errors handling
            df[column] = df[column].dropna().astype(dtype)
    except (ValueError, TypeError) as e:
        logging.error(f"Could not convert column '{column}' to type '{dtype}'.")
        logging.error(f"Conversion terminated unexpectedly with error: {e}.")
        logging.info(f"Column '{column}' appears to have dtype '{infer_dtype(column)}'")
        logging.info(
            f"Are you sure all values in {column} can be converted to type {dtype}?"
        )
        sys.exit(1)

    # Drop rows where conversion resulted in NaN (failed conversions)
    df = df.dropna(subset=[column])

    return df


def apply_column_types(
    df: pd.DataFrame, column_types: dict, error_on_invalid_type: bool = False
) -> pd.DataFrame:
    """
    Applies the given column-to-type mappings to the DataFrame.

    Parameters:
    df (pd.DataFrame): The DataFrame to modify.
    column_types (dict): A dictionary where keys are column names and values are target types.

    Returns:
    pd.DataFrame: The DataFrame with updated column types.
    """
    for col, dtype in column_types.items():
        if col in df.columns:
            df = safe_convert_column(df, col, dtype, error_on_invalid_type)
        else:
            logging.error(f"Column '{col}' not found in DataFrame.")
            sys.exit(1)

    return df


def set_up_logging(log_file) -> None:
    logfile = Path(".") / log_file
    fh = logging.FileHandler(logfile, mode="w")
    fh.setLevel(logging.DEBUG)
    sh = logging.StreamHandler()
    sh.setLevel(logging.INFO)
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[fh, sh],
    )


def main():
    args = parse_arguments()
    set_up_logging(args.logfile)

    # Parse input manifest into datatypes and filters
    parsed_manifest = parse_filter_manifest(args.filter_manifest)
    column_types = {}
    filters = []
    for column, filter_str, dtype in parsed_manifest.itertuples(index=False, name=None):
        column_types[column] = dtype
        filters.append(filter_str)

    # Read and clean the metadata TSV
    df = pd.read_csv(
        args.input,
        sep="\t",
        na_values=args.missing_values,
        skip_blank_lines=True,
        low_memory=False,
    )

    # Convert columns to appropriate types
    # (filters out rows where values do not convert)
    df = apply_column_types(
        df, column_types, error_on_invalid_type=args.error_on_invalid_type
    )

    # Filter using conditions
    df = apply_filters(df, filters)

    # Select output columns
    if args.select:
        df = df[args.select]

    # Output the filtered DataFrame
    if args.output:
        df.to_csv(args.output, sep="\t", index=False, header=(not args.remove_header))
        logging.info(f"Filtered data saved to {args.output}")
    else:
        # Print to stdout
        print(df)


if __name__ == "__main__":
    main()
