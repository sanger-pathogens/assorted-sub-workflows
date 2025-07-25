#!/usr/bin/env python3

import pandas as pd
import argparse
import logging
import sys
from pathlib import Path

pd.set_option("display.max_colwidth", None)
pd.set_option("display.max_rows", None)
pd.set_option("display.max_columns", None)

def setup_logging(log_file='qc_merge.log'):
    logging.basicConfig(
        filename=log_file,
        filemode='w',
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    logging.info("Logging initialized.")

def read_tsv(path, name):
    try:
        df = pd.read_csv(path, sep='\t')
        logging.info(f"{name} loaded successfully: {path}")
    except Exception as e:
        logging.error(f"Error reading {name} from {path}: {e}")
        sys.exit(1)
    return df

def process_gtdbtk(df):
    return df.rename(columns={
        'user_genome': 'genome_name_preqc',
        'classification': 'gtdbtk_taxonomy',
        'closest_genome_reference': 'closest_genome',
    })[[f'genome_name_preqc', 'gtdbtk_taxonomy', 'closest_genome']]


def process_checkm2(df, qc_stage):
    return df.rename(columns={
        'Name': f'genome_name_{qc_stage}',
        'Completeness': f'{qc_stage}_completeness',
        'Contamination': f'{qc_stage}_contamination',
        'Genome_Size': f'{qc_stage}_genome_size',
    })[[f'genome_name_{qc_stage}', f'{qc_stage}_completeness', f'{qc_stage}_contamination', f'{qc_stage}_genome_size']]

def process_gunc(df, qc_stage):
    return df.rename(columns={
        'genome': f'genome_name_{qc_stage}',
        'pass.GUNC': f'{qc_stage}_gunc_pass_or_fail',
        'n_contigs': f'{qc_stage}_contig_count'
    })[[f'genome_name_{qc_stage}', f'{qc_stage}_gunc_pass_or_fail', f'{qc_stage}_contig_count']]

def process_postqc_genome_name(df):
    new_df = df.copy() 
    new_df["genome_name_preqc"] = df["genome_name_postqc"].str.extract(r"cleaned_(.*)_filtered_kept_contigs")
    return new_df

def enrich_fields(df):
    df['sample_or_strain_name'] = df['genome_name_preqc'].str.rsplit("_", expand=True, n=1).loc[:, 0]
    df['genome_status'] = df['genome_name_preqc'].apply(lambda x: "mag" if "MAG" in x.upper() else "isolate")
    return df

def merge_all(*dfs):
    df, other_dfs = dfs[0], dfs[1:]
    join_keys = ['genome_name_preqc']  # specify your join keys
    for other_df in other_dfs:
        overlapping_cols = set(df.columns) & set(other_df.columns)
        cols_to_drop = overlapping_cols - set(join_keys)
        # print(other_df)
        other_df_clean = other_df.drop(columns=cols_to_drop)
        # print(other_df_clean)
        df = df.merge(other_df_clean, on=join_keys, how='left')
    # print(df)
    df = enrich_fields(df)
    # print(df)
    return df

def parse_args():
    parser = argparse.ArgumentParser(description="Merge GUNC, CheckM2 and post-QC metadata.")
    parser.add_argument('--pre_qc_checkm2', required=True, help='Path to pre-QC CheckM2 TSV')
    parser.add_argument('--pre_qc_gunc', required=True, help='Path to pre-QC GUNC TSV')
    parser.add_argument('--post_qc_checkm2', required=True, help='Path to post-QC CheckM2 TSV')
    parser.add_argument('--post_qc_gunc', required=True, help='Path to post-QC GUNC TSV')
    parser.add_argument('--gtdbtk', required=True, help='Path to GTDBTK TSV')
    parser.add_argument('--output', required=True, help='Output TSV path')
    return parser.parse_args()

def main():
    setup_logging()
    args = parse_args()

    gtdbtk = process_gtdbtk(read_tsv(args.gtdbtk, "GTDBTK classification"))
    pre_checkm2 = process_checkm2(read_tsv(args.pre_qc_checkm2, "Pre-QC CheckM2"), "preqc")
    pre_gunc = process_gunc(read_tsv(args.pre_qc_gunc, "Pre-QC GUNC"), "preqc")
    post_checkm2 = process_postqc_genome_name(process_checkm2(read_tsv(args.post_qc_checkm2, "Post-QC CheckM2"), "postqc"))
    post_gunc = process_postqc_genome_name(process_gunc(read_tsv(args.post_qc_gunc, "Post-QC GUNC"), "postqc"))

    merged = merge_all(pre_checkm2, pre_gunc, post_checkm2, post_gunc, gtdbtk)

    output_cols = [
        'genome_name_preqc',
        'genome_name_postqc',
        'sample_or_strain_name',
        'genome_status',
        'gtdbtk_taxonomy',
        'closest_genome',
        'preqc_genome_size',
        'preqc_contig_count',
        'preqc_completeness',
        'preqc_contamination',
        'preqc_gunc_pass_or_fail',
        'postqc_genome_size',
        'postqc_contig_count',
        'postqc_completeness',
        'postqc_contamination',
        'postqc_gunc_pass_or_fail',
        'mdm_cleaner_status',
        'final_qc_status',
    ]
    # print(merged)
    for col in output_cols:
        if col not in merged.columns:
            merged[col] = 'NA'
    # print(merged)

    merged[output_cols].to_csv(args.output, sep='\t', index=False)
    logging.info(f"Merged QC report written to: {args.output}")

if __name__ == '__main__':
    main()
