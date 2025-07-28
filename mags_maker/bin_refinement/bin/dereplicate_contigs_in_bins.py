#!/usr/bin/env python3

import sys
import os
from collections import defaultdict
import logging

'''
Usage:
    ./dereplicate_contigs_in_bins.py checkm2_report.txt input_dir output_dir ["remove_ambiguous_contigs"]

Arguments:
    checkm2_report.txt  Tab-seperated report from CheckM2 containing Completeness and Contamination scores
    input_dir                   Directory with input reads as FASTAs (bin file?)
    output_dir                  Directory to output dereplicated bins
    remove_ambiguous_contigs    Optional (string), remove contigs shared by multiple bins (strict bin purity)
'''

# Set up logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)

# Score bins using completeness and contamination scores from CheckM2 report
logging.info("Loading in bin completeness and contamination scores...")
bin_scores = {}

required_columns = {"Completeness", "Contamination"}

try:
    with open(sys.argv[1]) as bin_file:
        # Check the header on the report is as expected
        first_line = next(bin_file, None)
        if first_line is None:
            raise ValueError(f"Error: File {sys.argv[1]} is empty.")
        
        header_columns = set(first_line.strip().split("\t"))
        missing_columns = required_columns - header_columns
        if missing_columns:
            raise ValueError(f"Header in {sys.argv[1]} is missing required columns: {missing_columns}")

        if not required_columns.issubset(header_columns):
            raise ValueError(f"Header of {sys.argv[1]} is missing the required columns: {required_columns}")
        
        # Then process each non-header line to determine scores per bin
        for line in bin_file:
            cut = line.strip().split("\t")   #TODO: Move to pandas? Currently reliant on consistent checkm2 reports
            score = float(cut[1]) - 5 * float(cut[2]) + 0.000_000_0001 * float(cut[5])
            bin_scores[cut[0]] = score

except FileNotFoundError as e:
    logging.error(f"Error: The file {sys.argv[1]} was not found.\n With error: {e}")
    sys.exit(1)

except IOError as e:
    logging.error(f"Error: An I/O error occurred while reading {sys.argv[1]}.\n With error: {e}")
    sys.exit(1)

except Exception as e:
    logging.exception(f"An unexpected error occurred: {e}")
    sys.exit(1)

# Assign contigs to the highest scoring bin they appear in, unless remove_ambiguous_contigs option is given
logging.info("Loading in contigs in each bin...")
contig_mapping = defaultdict(str)
try:
    for bin_file in os.listdir(sys.argv[2]):
        bin_name = ".".join(bin_file.split("/")[-1].split(".")[:-1])
        with open(os.path.join(sys.argv[2], bin_file)) as f:
            for line in f:
                if not line.startswith(">"):    #TODO: move to Biopython for safer and easier fasta handling?
                    continue
                contig = line[1:].strip()
                if contig not in contig_mapping:
                    contig_mapping[contig] = bin_name
                else:
                    if len(sys.argv) > 4 and sys.argv[4] == "remove_ambiguous_contigs":
                        contig_mapping[contig] = None
                    elif bin_scores.get(bin_name, 0) > bin_scores.get(contig_mapping[contig], 0):
                        contig_mapping[contig] = bin_name

except IOError as e:
    logging.error(f"An IO error has occured with the file {bin_file}.\nexiting with error: {e}")
    sys.exit(1)
except FileNotFoundError as e:
    logging.error(f"The program was unable to find the file {bin_file}.\nexiting with error: {e}")
    sys.exit(1)
except Exception as e:
    logging.exception(f"An unknown error has occured: {e}")
    sys.exit(1)

# Make a new dereplicated version of each bin file based on the final assignments
logging.info("Making a new dereplicated version of each bin...")
os.makedirs(sys.argv[3], exist_ok=True)

for bin_file in os.listdir(sys.argv[2]):
    bin_name = ".".join(bin_file.split("/")[-1].split(".")[:-1])
    input_path = os.path.join(sys.argv[2], bin_file)
    output_path = os.path.join(sys.argv[3], bin_file)
    
    with open(input_path) as infile, open(output_path, 'w') as out:
        at_least_one = False
        store = False
        for line in infile:
            if line.startswith(">"):
                contig = line[1:].strip()
                if contig_mapping.get(contig) == bin_name:
                    at_least_one = True
                    store = True
                else:
                    store = False
            if store:
                out.write(line)

    # Remove any bins that are empty by the end
    if not at_least_one:
        os.remove(output_path)
