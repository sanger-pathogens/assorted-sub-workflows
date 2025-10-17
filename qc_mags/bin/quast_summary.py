#!/usr/bin/env python3

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd


def setup_logging(log_file: str = "quast_summary.log"):
    logging.basicConfig(
        level=logging.INFO,
        handlers=[logging.StreamHandler(), logging.FileHandler(log_file, mode="w")],
        format="%(asctime)s - %(levelname)s - %(message)s",
    )
    logging.info("Logging initialized.")


def read_quast_tsv(file_path: Path):
    try:
        quast_df = pd.read_csv(file_path, sep='\t')
        logging.info(f'File read successfully -> {file_path}')
    except IOError as e:
        logging.error(f'Error reading file: {file_path} -> {e}')
        sys.exit(1)
    except Exception as e:
        logging.error(f'An unknown exception occurred: {file_path} -> {e}')
        sys.exit(1)
    return quast_df



def calculate_small_contigs(df):
    try:
        small_contig_df = df[['Assembly', 'N50', 'Total length', '# contigs (>= 0 bp)', '# contigs (>= 1000 bp)']]
    except Exception as e:
        logging.error(f'Unable to parse headers -> {e}')
        sys.exit(1)
    
    try:
        small_contig_df['small contigs'] = (small_contig_df['# contigs (>= 0 bp)'] - small_contig_df['# contigs (>= 1000 bp)'])
        small_contig_df['Proportion Contigs <= 1kbp'] = (small_contig_df['small contigs'] / small_contig_df['# contigs (>= 0 bp)'])
        logging.info('Computed proportion of small contigs')
    except Exception as e:
        logging.error(f'Unable to compute small contigs -> {e}')
        sys.exit(1)
    return small_contig_df


def write_output_tsv(small_contig_df, out_path):
    try:
        summary_df = small_contig_df[['Assembly', 'N50', 'Total length', 'Proportion Contigs <= 1kbp']]
        summary_df.to_csv(out_path, sep='\t', index=False)
        logging.info(f'Successfully created summary file -> {out_path}')
    except Exception as e:
        logging.error(f'Failed to create summary file -> {e}')
        sys.exit(1)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Summary quality of contigs from quast stats"
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help='Path to transposed_report.tsv quast output'
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

    small_contigs = calculate_small_contigs(read_quast_tsv(args.input))
    write_output_tsv(small_contigs, args.output)


if __name__ == "__main__":
    main()
