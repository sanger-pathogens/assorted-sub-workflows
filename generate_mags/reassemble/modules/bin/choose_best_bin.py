#!/usr/bin/env python3

import argparse
import csv
import logging
import os
import shutil
from pathlib import Path

def parse_arguments():

    def restricted_float(x):
        try:
            x = float(x)
        except ValueError:
            raise argparse.ArgumentTypeError(f"{x} not a floating-point literal")
        
        if x < 0.0 or x > 100.0:
            raise argparse.ArgumentTypeError(f"{x} not in range [0.0, 100.0]")
        return x

    parser = argparse.ArgumentParser(
        description="Select the best bins based on completeness and contamination thresholds."
    )
    parser.add_argument("checkm_tsv", help="CheckM output TSV file")
    parser.add_argument("bin_dir", help="Directory containing bin FASTA files")
    parser.add_argument("output_dir", help="Directory to save the best bins")
    parser.add_argument("--min-completeness", type=restricted_float, default=50.0,
                        help="Minimum completeness threshold (default: 50.0)")
    parser.add_argument("--max-contamination", type=restricted_float, default=5.0,
                        help="Maximum contamination threshold (default: 5.0)")
    parser.add_argument("--log", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help="Logging level (default: INFO)")
    return parser.parse_args()

def setup_logging(level):
    logging.basicConfig(format="%(levelname)s: %(message)s", level=getattr(logging, level))

def score_bin(completeness, contamination):
    return completeness + 5 * (100 - contamination)

def select_best_bins(tsv_file, min_compl, max_contam):
    best_bins = {}
    with open(tsv_file, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            completeness = float(row['completeness'])
            contamination = float(row['contamination'])

            if completeness < min_compl or contamination > max_contam:
                continue

            bin_full = row['bin']
            bin_parts = bin_full.split('_')
            style = bin_parts[-1]
            bin_base = '_'.join(bin_parts[:-1])

            score = score_bin(completeness, contamination)
            n50 = int(row['N50'])

            if bin_base not in best_bins or (
                score > best_bins[bin_base][1] or 
                (score == best_bins[bin_base][1] and n50 > best_bins[bin_base][2])
            ):
                best_bins[bin_base] = (style, score, n50)

    return best_bins

def copy_best_bins(best_bins, bin_dir, output_dir):
    ''' This function copy the best bins to a destination folder '''
    os.makedirs(output_dir, exist_ok=True)
    copied = 0
    for bin_base, (style, _, _) in best_bins.items():
        bin_filename = f"{bin_base}_{style}.fasta"
        bin_path = Path(bin_dir) / bin_filename
        if not bin_path.exists():
            logging.warning(f"Bin file not found: {bin_filename}")
            continue
        dest_path = Path(output_dir) / bin_path.name
        shutil.copy(bin_path, dest_path)
        logging.info(f"Copied: {bin_path} -> {dest_path}")
        copied += 1
    logging.info(f"Total bins copied: {copied}")

    if len(bin_base.keys())!= copied:
        logging.error(f"Mismatch between number of bins and number of bins copied.\nExpected to copy {len(bin_base.keys())} bins but instead copied {copied}")

def write_filtered_tsv(original_tsv, best_bins, output_tsv):
    with open(original_tsv, newline='') as infile, open(output_tsv, 'w', newline='') as outfile:
        reader = csv.DictReader(infile, delimiter='\t')
        fieldnames = reader.fieldnames
        writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()

        for row in reader:
            bin_full = row['bin']
            bin_parts = bin_full.split('_')
            style = bin_parts[-1]
            bin_base = '_'.join(bin_parts[:-1])

            if bin_base in best_bins and best_bins[bin_base][0] == style:
                writer.writerow(row)
    logging.info(f"Filtered TSV written to: {output_tsv}")

def main():
    args = parse_arguments()
    setup_logging(args.log)

    logging.info("Selecting best bins...")
    best_bins = select_best_bins(args.checkm_tsv, args.min_completeness, args.max_contamination)

    logging.info(f"Found {len(best_bins)} best bins. Copying to {args.output_dir}...")
    copy_best_bins(best_bins, args.bin_dir, args.output_dir)

    filtered_tsv = Path(f"{args.output_dir}_summary.tsv")
    write_filtered_tsv(args.checkm_tsv, best_bins, filtered_tsv)

if __name__ == "__main__":
    main()
