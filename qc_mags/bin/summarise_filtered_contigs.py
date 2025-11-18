#!/usr/bin/env python3

"""
Generate a quick summary of number and proportion of contigs filtered out by filtering stages.
Will need to be adapted if further contig filtering steps are included: only works for length 
filter (currently implemented by seqkit), and contamination filter (mdmcleaner).

Usage:
    summarise_filtered_contigs.py
        --preqc_quast <quast_report1.tsv> \
        --postqc_quast <quast_report2.tsv> \
        --min_contig <int> \
        --output <output/path>
"""

import argparse
import pandas as pd
from pathlib import Path
import logging
import sys

def main():
    args = parse_args()
    setup_logging(args.debug)

    preqc_data = read_tsv(args.pre_qc_quast, "preqc QUAST report")
    postqc_data = read_tsv(args.post_qc_quast,  "postqc QUAST report")
    postqc_data = prep_postqc_join_column(preqc_data, postqc_data)

    data = join_quast(
        preqc_df = preqc_data, 
        postqc_df = postqc_data, 
        min_contig = args.min_contig
        )

    enriched_data = generate_output_columns(min_contig=args.min_contig, merged_df=data)
    
    #write_tsv(enriched_data, args.output, "filtering_summary")
    write_tsv(enriched_data, args.output)

def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate a quick summary of number and proportion of contigs filtered out by which stage.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--debug", 
        type=bool, 
        required=False, 
        default=False,
        help="Set to true to log debug-level information"
    )
    parser.add_argument(
        "--pre_qc_quast",
        type=Path,
        required=True,
        help="QUAST report (transposed TSV format) run on unfiltered assemblies."
    )
    parser.add_argument(
        "--post_qc_quast",
        type=Path,
        required=True,
        help="QUAST report (transposed TSV format) run on assemblies after filtering contigs."
    )
    parser.add_argument(
        "--min_contig",
        type=int,
        required=True,
        help="Length (in bp) below which contigs have been filtered out."
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output TSV path."
    )
    return parser.parse_args()

def setup_logging(debug: bool):
    log_filename = "summarise_filtered_contigs.log"
    
    if debug:
        logging.basicConfig(
        level=logging.DEBUG,
        handlers=[logging.StreamHandler(), logging.FileHandler(log_filename, mode="w")],
        format="%(asctime)s - %(levelname)s - %(message)s",
    )
    else:
        logging.basicConfig(
        level=logging.INFO,
        handlers=[logging.StreamHandler(), logging.FileHandler(log_filename, mode="w")],
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
        
def check_length_filter(length_fail_col:str, postqc_df:pd.DataFrame):
    # Warn if post qc bins contain contigs below min_contig
    failing_mask = postqc_df[length_fail_col] > 0    # negative value equal number small contigs that also fail mdmcleaner
    if failing_mask.any():
        failing_assemblies = postqc_df.loc[failing_mask, 'Assembly'].tolist()
        logging.warning(f"Contigs below the supplied minimum contig length are found by QUAST in postqc bins: "
                        f"{failing_assemblies}. Check seqkit is functioning correctly.")

def prep_postqc_join_column(preqc_df: pd.DataFrame, postqc_df: pd.DataFrame):
    # 'Assembly' name for pre and postqc need to be identical to join on this col downstream

    postqc_df = postqc_df.rename(columns={'Assembly':'postqc_Assembly'})

    # New Assembly column with values from matched preqc names
    for i, row in postqc_df.iterrows():
        postqc_name = row['postqc_Assembly']
        found_match = False

        for preqc_name in preqc_df['Assembly']:
            if preqc_name in postqc_name:
                postqc_df.at[i, 'Assembly'] = preqc_name
                found_match = True
                break
        
        if not found_match:
            logging.warning("Cannot match pre- and post-qc QUAST report based on values in Assembly columns. \
                            Expected preqc name to be a substring of postqc name.")

    # Log matched assembly names pre and post qc for debugging:
    for _, row in postqc_df.dropna(subset=['Assembly']).iterrows():
        logging.debug(f"Matched preqc Assembly: {row['Assembly']}  -->  postqc Assembly: {row['postqc_Assembly']}")
    
    return postqc_df

def join_quast(preqc_df:pd.DataFrame, postqc_df:pd.DataFrame, min_contig:int) -> pd.DataFrame:
    # Add columns for small contigs]
    length_pass_col = f"# contigs (>= {min_contig} bp)"
    if length_pass_col not in preqc_df.columns or length_pass_col not in postqc_df.columns:
        logging.error(f"Expected column {length_pass_col} not found in DataFrame columns")

    length_fail_col = f"# contigs (< {min_contig} bp)"
    preqc_enriched_df = calc_small_contigs(length_pass_col, length_fail_col, preqc_df)
    postqc_enriched_df = calc_small_contigs(length_pass_col, length_fail_col, postqc_df)

    # Double check that filtering has taken place, warning if not
    check_length_filter(length_fail_col, postqc_enriched_df)

    # Join on Assembly ID, bringing only necessary columns into the merged df
    required_cols = ['Assembly', '# contigs', f'# contigs (< {min_contig} bp)']
    merged_df = pd.merge(
        preqc_enriched_df[required_cols],
        postqc_enriched_df[required_cols + ['postqc_Assembly']],
        on=['Assembly'],
        how='outer',    # Will put NaNs in any missing fields (i.e. postqc name in preqc data)
        suffixes=('_preqc', '_postqc')    # avoids identical col clashes
    )

    # TODO: debugging step printing any NaN containing rows?

    return merged_df

def calc_small_contigs(length_pass_col:str, length_fail_col:str, df:pd.DataFrame) -> pd.DataFrame:
    # Add column for number of small contigs
    enriched_df = df.copy()
    enriched_df[length_fail_col] = df['# contigs'] - df[length_pass_col]
    return enriched_df

def generate_output_columns(merged_df:pd.DataFrame, min_contig:int) -> pd.DataFrame:
    # Add cols for number and percent of contigs filtered out due to length, failing mdmcleaner and overall
    merged_df['total_contigs_filtered'] = merged_df['# contigs_preqc'] - merged_df['# contigs_postqc']
    merged_df['%_contigs_filtered'] = (merged_df['total_contigs_filtered'] * 100 / merged_df['# contigs_preqc']).round(2)


    merged_df['num_small_contigs'] = merged_df[f'# contigs (< {min_contig} bp)_preqc']
    merged_df['%_small_contigs'] = (merged_df['num_small_contigs'] * 100 / merged_df['# contigs_preqc']).round(2)

    merged_df['num_mdmcleaner_failed'] = merged_df['total_contigs_filtered'] - merged_df['num_small_contigs']
    merged_df['%_mdmcleaner_failed'] = (merged_df['num_mdmcleaner_failed'] * 100 / merged_df['# contigs_preqc']).round(2)

    final_df = merged_df[
        [
            'Assembly', 
            'postqc_Assembly',
            'total_contigs_filtered',
            '%_contigs_filtered',
            'num_small_contigs',
            '%_small_contigs',
            'num_mdmcleaner_failed',
            '%_mdmcleaner_failed'
        ]
    ].copy()

    # Rename assembly name columns to match report.py convention 
    final_df.rename(columns={'Assembly':'preqc_genome_name', 'postqc_Assembly':'postqc_genome_name'}, inplace=True)

    return final_df

# def write_tsv(final_df: pd.DataFrame, output_path:Path, filename:str):
#     output_file = output_path / f"{filename}.tsv"
#     final_df.to_csv(output_file, sep="\t", index=False)
#     logging.info(f"Contig {filename}.tsv written to: {output_path}")

def write_tsv(df: pd.DataFrame, output_path:Path):
    df.to_csv(output_path, sep="\t", index=False)
    logging.info(f"Contig filtering summary written to: {output_path}")

if __name__ == "__main__":
    main()