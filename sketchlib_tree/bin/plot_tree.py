#!/usr/bin/env python3

import argparse

from ete3 import Tree, TreeStyle


def plot_tree(newick_file, output_file):
    t = Tree(newick_file)
    ts = TreeStyle()
    t.render(output_file, w=183, units="mm", tree_style=ts)


def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Save a phylogenetic tree to a PNG file from a Newick file.")
    parser.add_argument("newick_file", type=str, help="Path to the Newick file.")
    parser.add_argument("output_file", type=str, help="Path to the output PNG file.")
    args = parser.parse_args()

    # Plot the tree
    plot_tree(args.newick_file, args.output_file)


if __name__ == "__main__":
    main()
