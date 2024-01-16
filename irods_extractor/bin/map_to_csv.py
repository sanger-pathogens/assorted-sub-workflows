#!/usr/bin/env python3

import argparse
import pandas as pd

#This file needs to be moved to the /bin/ folder inside the project directory intending to use this script.

def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="tool to write csv from nextflow map",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input_map_list",
        help="A collection of maps from nextflow to be written to csv",
    )
    return parser.parse_args()

def element_in_list_to_dict(input_list: list):
    output_list = []
    for i in input_list:
        metadata = {}
        pairs = i.split(', ')
        for pair in pairs:
            full_data= pair.split('=')
            key = full_data[0]
            value = full_data[1:]
            metadata[key] = "".join(value)
        output_list.append(metadata)
    return output_list
        

def dataframe_from_input_list(input_list: list):
   df = pd.DataFrame(input_list)
   df = df.set_index('ID')
   return df

args = parse_arguments()

with open(args.input_map_list) as file:
    input_data = file.readlines()

data_list = element_in_list_to_dict(input_data)

df = dataframe_from_input_list(data_list)

df.to_csv("metadata.csv")
