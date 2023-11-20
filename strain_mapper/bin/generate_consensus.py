#!/usr/bin/env python3

import argparse
import csv
import logging
from pathlib import Path
from typing import TextIO

DEFAULTSEQCHAR = 'N'

def get_chrom_id_and_size(ref_index: Path):
    with open(ref_index, "r") as f:
        seq_id, size, _ = f.readline().split(maxsplit=2)
    return seq_id, int(size)


def initialise_seq(chrom_size: int):
    return [DEFAULTSEQCHAR] * chrom_size


def parse_lines(fh: TextIO):
    for line in fh:
        if not line.startswith("##"):
            header = line.lstrip("#").split()
            break
    return csv.DictReader(fh, delimiter="\t", fieldnames=header)


def parse_position(pos: str):
    try:
        pos = int(pos)  # reference nucleotide position converted to int
    except ValueError:
        logging.warning(
            f"Expected integer value for nucleotide position but got: '{pos}'"
        )
        return None
    else:
        return pos


def parse_quality(qual: str):
    try:
        qual = float(qual)
    except ValueError:
        if qual == ".":
            return None
        else:
            raise
    else:
        return qual


def is_acceptable_quality(qual: float, pos: int, threshold: int=10, alt_quality: dict=None):
    if alt_quality is None:
         alt_quality = {}
    if qual is None:
        return False
    if qual > threshold or (pos in alt_quality and qual > alt_quality[pos]):
        return True
    return False


def get_nt_to_add(ref: str, alt: str):
    if len(alt) > 1:
        # if mutiple bases called at a position
        return DEFAULTSEQCHAR
    if alt == ".":
        # if the mapped strain is the same as the query, then it is reported as a '.'
        return ref
    # for SNPs
    return alt


def get_seq(
    vcf: Path, ref_index: Path, qual_threshold: float=10
):
    chrom_id, chrom_size = get_chrom_id_and_size(ref_index)
    seq = initialise_seq(chrom_size)
    with open(vcf, "r", newline="") as f:
        parsed_lines = parse_lines(f)
        for line in parsed_lines:
            seq_id = line["CHROM"]
            ref = line["REF"]  # reference base
            alt = line["ALT"]  # alternative base
            pos = parse_position(line["POS"])
            qual = parse_quality(line["QUAL"])
            if seq_id == chrom_id:
                if qual is None:
                    logging.warning(
                        f"The following line had an unexpected quality value: {line}"
                    )
                if is_acceptable_quality(
                    qual, pos, qual_threshold
                ):
                    seq[pos - 1] = get_nt_to_add(ref, alt)
    assert len(seq) == chrom_size
    return "".join(seq)


def write_seq(seq: str, seq_id: str, output_fasta: str="out.fa"):
    with open(output_fasta, "w") as f:
        seq_header = f">{seq_id}"
        f.write(f"{seq_header}\n{seq}\n")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Script to align map reads back to reference"
    )
    parser.add_argument(
        "--sample",
        "-s",
        help="Sample ID of sample whose sequence will be aligned to the reference",
    )
    parser.add_argument(
        "--vcf", "-v", help="VCF (generated by `bcftools call`) for the sample"
    )
    parser.add_argument(
        "--ref-index",
        "-i",
        help="Reference index (.fai) for the strain to which the sample will be mapped",
    )
    parser.add_argument(
        "--output",
        "-o",
        help="Output filename for aligned reads (fasta format)",
        default="out.fa",
    )
    return parser.parse_args()


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    args = parse_args()

    seq = get_seq(args.vcf, args.ref_index)
    write_seq(seq, args.sample, args.output)