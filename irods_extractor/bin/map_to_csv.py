#!/usr/bin/env python3

import argparse
import json
import pandas as pd

#This file needs to be moved to the /bin/ folder inside the project directory intending to use this script.

def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Tool to write csv from nextflow maps serialised as JSON Lines",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input_map_list",
        help="A collection of maps from nextflow, one JSON object per line",
    )

    parser.add_argument("--output",
        default="metadata.csv",
        help="file to output data to",
    )
    return parser.parse_args()

def records_from_jsonl(input_lines: list, source: str):
    """Parse one map per line.

    The maps are serialised with JsonOutput.toJson rather than Map.toString so
    that values containing ',', '=' or a newline survive the round trip.
    """
    records = []
    for line_number, line in enumerate(input_lines, start=1):
        if not line.strip():
            continue
        try:
            records.append(json.loads(line))
        except json.JSONDecodeError as error:
            raise ValueError(
                f"{source} line {line_number} is not a JSON object: {error}"
            ) from error
    return records


def dataframe_from_input_list(input_list: list):
   df = pd.DataFrame(input_list)
   df = df.set_index('ID')
   return df

args = parse_arguments()

with open(args.input_map_list) as file:
    input_data = file.readlines()

data_list = records_from_jsonl(input_data, args.input_map_list)

df = dataframe_from_input_list(data_list)

df.to_csv(args.output)
