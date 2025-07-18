#!/usr/bin/env python3
import sys
import os
from collections import defaultdict
import logging

# Usage: ./script.py bins.stats binsFolder outFolder

# Load in bin completion and contamination scores
print("Loading in bin completion and contamination scores...")
bin_scores = {}
try:
    with open(sys.argv[1]) as bin_file:
        for line in bin_file:
            if "Completeness" in line: # Skip the header
                continue
            cut = line.strip().split("\t")
            score = float(cut[1]) - 5 * float(cut[2]) + 0.000_000_0001 * float(cut[5])
            bin_scores[cut[0]] = score

except FileNotFoundError as e:
    print(f"Error: The file {sys.argv[1]} was not found.\n With error: {e}")
    sys.exit(1)

except IOError as e:
    print(f"Error: An I/O error occurred while reading {sys.argv[1]}.\n With error: {e}")
    sys.exit(1)

except Exception as e:
    print(f"An unexpected error occurred: {e}")
    sys.exit(1)

# Load in contigs in each bin
print("Loading in contigs in each bin...")
contig_mapping = defaultdict(str)
try:
    for bin_file in os.listdir(sys.argv[2]):
        bin_name = ".".join(bin_file.split("/")[-1].split(".")[:-1])
        with open(os.path.join(sys.argv[2], bin_file)) as f:
            for line in f:
                if not line.startswith(">"):
                    continue
                contig = line[1:].strip()
                if contig not in contig_mapping:
                    contig_mapping[contig] = bin_name
                else:
                    if len(sys.argv) > 4 and sys.argv[4] == "remove":
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
    logging.error(f"An unknown error has occured: {e}")
    sys.exit(1)

# Go over the bin files again and make a new dereplicated version of each bin file
print("Making a new dereplicated version of each bin file")
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
    
    if not at_least_one:
        os.remove(output_path)
