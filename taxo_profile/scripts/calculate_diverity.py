#!/usr/bin/env python3

import os
from pathlib import Path
import argparse
from skbio.diversity import alpha_diversity
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description="A Basic tool to calculate diversity of a metagenomic sample using abundance estimates derived from sylph output")
    parser.add_argument("-i", "--sylph-summary-input", nargs="+", required=True, help="A sylph summary file derived from the taxo_profile subworkflow")
    parser.add_argument("-a", "--taxonomic-abundance-threshold", type=float, required=False, help="Only calculate diversity based on the samples with taxonomic_abundance >= taxonomic-abundance-threshold")
    parser.add_argument("-o", "--outdir", type=Path, default=Path.cwd(), help="")
    return parser.parse_args() 

def calculate_alpha_sylph(abund_matrix):
    """
    Calculates alpha diversity from sylph abundanc estimate
    """
    counts = abund_matrix.values
    sample_ids = abund_matrix.index.tolist()

    # Shannon diversity
    shannon = alpha_diversity(metric="shannon", counts=counts, ids=sample_ids)

    # Simpson's index
    simpson = alpha_diversity(metric="simpson", counts=counts, ids=sample_ids)

    # Observed richness (number of taxa with abundance > 0)
    richness = alpha_diversity(metric="observed_otus", counts=counts, ids=sample_ids)

    # Pielou's evenness
    pielou = alpha_diversity(metric="pielou_e", counts=counts, ids=sample_ids)

    berger_parker = alpha_diversity(metric="berger_parker_d", counts=counts, ids=sample_ids)

    results = pd.DataFrame({
        "shannon": shannon,
        "simpson": simpson,
        "richness": richness,
        "pielou_evenness": pielou, 
        "berger_parker_dominance":berger_parker
    })

    results.index.name = "ID"
    return results

    

def load_data(data_summary: list, abundance_threshold=None):
    abund_matrices = []
    for sylph in data_summary:
        df = pd.read_csv(sylph, sep="\t")
        if not abundance_threshold:
            pass
        else:
            df = df[df["Taxonomic_abundance"] >= abundance_threshold]
        abund_matrix = df.pivot_table(
            index="Sample_file",
            columns="Genome_file",
            values="Taxonomic_abundance",
            fill_value=0
        )
        abund_matrices.append(abund_matrix)
    combined_matrix = pd.concat(abund_matrices, axis=0).fillna(0)

    return combined_matrix 

def main():
    args = parse_args()

    if not args.taxonomic_abundance_threshold:
        abund_matrix = load_data(args.sylph_summary_input)
    else: 
        abund_matrix = load_data(args.sylph_summary_input, abundance_threshold=args.taxonomic_abundance_threshold)
    
    results = calculate_alpha_sylph(abund_matrix)
    results.to_csv(args.outdir / "diversity_estimates.tsv", sep='\t')

if __name__ == "__main__":
    main()