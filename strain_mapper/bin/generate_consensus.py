#!/usr/bin/env python3

import argparse
import csv
import logging
from pathlib import Path
from typing import TextIO
from collections import OrderedDict, Counter

def log_filter_info(filter_counter: Counter, used_counter: Counter):
    """
    Logs information about FILTER values encountered and used during VCF parsing,
    including the percentage of ambiguous variants.

    Parameters:
    filter_counter (Counter): A Counter object tallying all FILTER field values
        encountered in the VCF file (e.g., PASS, ., lowQual, het).
    used_counter (Counter): A Counter object tallying only the FILTER values that
        were actually used to incorporate variants into the final consensus sequence
        (typically only PASS and/or .).
    total_variants (int): The total number of variants processed from the VCF file.

    Behavior:
    - Logs all unique FILTER values found in the VCF and their frequencies.
    - Warns if both 'PASS' and '.' are present in the VCF, which can imply mixed
      filtering conventions.
    - Reports basic stats on consensus i.e how many ambigious bases
      
    Notes:
    - Assumes that only variants with FILTER == 'PASS' or '.' are used.
    """
    all_filters = sorted(filter_counter.keys())
    logging.info(f"Found FILTER values in VCF: {', '.join(all_filters)}")
    
    if "PASS" in filter_counter and "." in filter_counter:
        logging.warning("VCF contains both 'PASS' and '.' in the FILTER column. "
                        "These are treated as equivalent, but consider changing.")
    
    # this will be nice to show if you have just 1 random dot somewhere
    if "PASS" in filter_counter:
        logging.info(f"Variants used for consensus (FILTER == 'PASS'): {used_counter['PASS']}")
    if "." in filter_counter:
        logging.info(f"Variants used for consensus (FILTER == '.'): {used_counter['.']}")
    
    ambiguous = sum(used_counter[f] for f in used_counter if f not in {"PASS", "."})
    if ambiguous > 0:
        ambiguous_percentage = (ambiguous / total_variants) * 100
        logging.info(f"Variants with non-passing FILTERs (left ambiguous): {ambiguous} "
                     f"({ambiguous_percentage:.2f}% of total variants)")


def get_chrom_id_and_size(ref_index: Path) -> dict[str, int]:
    """
    Extracts chromosome IDs and their corresponding sizes from a reference index file.

    Parameters:
    ref_index (Path): Path to the reference index file generated by samtools (or similar).

    Returns:
    dict[str, int]: A dictionary containing chromosome IDs as keys and their corresponding sizes as values.

    This function only consideres the two first elements in the split which are name and length

    Example:
    If the reference index file contains:
    ```
    NZ_CP012480_1_1	2074179	
    NZ_CP012742_1_2	4944	
    ```
    The function will return:
    ```
    {'NZ_CP012480_1_1': 2074179, 'NZ_CP012742_1_2':	4944}
    """
    chrom_id_size = OrderedDict()
    with open(ref_index, "r") as f:
        for line in f:
            seq_id, size, _ = line.split(maxsplit=2)
            chrom_id_size[seq_id] = int(size)
    return chrom_id_size


def initialise_seq(chrom_size: int, default_seq_character) -> list[str]:
    """
    Initializes a list representing a sequence with default characters.

    Parameters:
    chrom_size (int): The size of the chromosome sequence.

    Example:
    Tasked with a genome size of 2_000_000 and a default character of N
    this function will return a list of N 2_000_000 long
    """
    return [default_seq_character] * chrom_size


def parse_lines(fh: TextIO) -> csv.DictReader:
    """
    Parses lines from a Variant Call Format (VCF) file handle and returns a CSV DictReader.

    Parameters:
    fh (TextIO): A file handle to read lines from.

    Returns:
    csv.DictReader: A DictReader object for the parsed CSV data.

    This function reads lines from the provided VCF file handle. It ignores lines starting
    with '##'.

    It extracts the header from the first non-comment line of the VCF file, which
    typically contains column names such as 'POS', 'QUAL', etc.

    It then creates a DictReader object using the remaining lines in the file handle, treating
    them as tab-separated data, and returns it.

    Example:
    Suppose the VCF file handle contains:
    ```
    ##fileformat=VCFv4.3
    ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
    #CHROM  POS ID  REF ALT QUAL    FILTER  INFO
    chr1    1000    .   T   A   10.5    PASS    AF=0.3
    chr2    2000    .   G   C   20.7    PASS    AF=0.1
    ```
    The function will return a DictReader object configured to read CSV data with the
    specified header ('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'):
    ```
    <csv.DictReader object at 0x7f84b31aa310>
    ```
    """
    for line in fh:
        if not line.startswith("##"):
            header = line.lstrip("#").split()
            break
    return csv.DictReader(fh, delimiter="\t", fieldnames=header)


def parse_position(pos: str) -> int:
    """
    Parses a string representing a nucleotide position and returns an integer.

    Parameters:
    pos (str): A string representing a nucleotide position.

    Returns:
    An integer representing the parsed nucleotide position,
    or None if the input string cannot be converted to an integer.

    Example:
    If pos is '123', the function will return:
    123

    If pos is 'abc', the function will log a warning and return:
    None
    """
    try:
        pos = int(pos)  # reference nucleotide position converted to int
    except ValueError:
        logging.warning(
            f"Expected integer value for nucleotide position but got: '{pos}'"
        )
        return None
    else:
        return pos


def get_nt_to_add(ref: str, alt: str, default_seq_character):
    """
    Determines the nucleotide to add to the sequence based on the reference and alternative alleles.

    Parameters:
    ref (str): The reference allele.
    alt (str): The alternative allele.
    default_seq_character (str): The default nucleotide character to use. (N)

    Returns:
    str: The nucleotide character to add to the sequence.

    This function determines the nucleotide to add to the sequence based on the reference allele 'ref',
    the alternative allele 'alt', and the default nucleotide character 'default_seq_character'.
    If the length of 'alt' is greater than 1, indicating multiple bases called at a position,
    it returns the default nucleotide character.
    If 'alt' is '.', indicating that the mapped strain is the same as the query, it returns 'ref'.
    Otherwise, for single nucleotide polymorphisms (SNPs), it returns 'alt'.

    Example:
    If ref is 'A', alt is 'T', and default_seq_character is 'N', the function will return:
    'T'

    If ref is 'C', alt is '.', and default_seq_character is 'N', the function will return:
    'C'

    If ref is 'G', alt is 'CG', and default_seq_character is 'N', the function will return:
    'N'
    """
    if len(alt) > 1:
        # if mutiple bases called at a position
        return default_seq_character
    if alt == ".":
        # if the mapped strain is the same as the query, then it is reported as a '.'
        return ref
    # for SNPs
    return alt


def get_seq(vcf: Path, ref_index: Path, default_seq_character: str) -> dict[str, str]:
    """
    Constructs sequences from VCF data, optionally filtering by FILTER status.

    Parameters:
    vcf (Path): Path to the input VCF file containing variant data.
    ref_index (Path): Path to the reference index file containing chromosome sizes.
    default_seq_character (str): Default character to initialize sequences with which 
        will appear at sites where no base has been called i.e. how to represent when neither a REF 
        nor a VARIANT call.

    Returns:
    Dict[str, str]: A dictionary where keys are chromosome IDs and values are the
        constructed sequences with variants incorporated.

    Notes:
    - variants are included if FILTER column contains "PASS" or "."
    - The function preserves VCF convention where empty FILTER ('.') indicates
      the variant has passed all filters

    Example:
    If vcf contains:
        chr1 100 A T PASS
        chr1 200 G C .
        chr1 300 C A FAIL
    Then:
    - get_seq(vcf, ref_index, 'N') would incorporate chr1:100 and chr1:200 into the sequence,
      but skip chr1:300 due to the failing filter.

    It is expected that positions are FILTERED consistently with either PASS or . being used,
    depending on the filtering approach used (soft vs. hard filtering).
    """
    chrom_id_size = get_chrom_id_and_size(ref_index)
    seq = {}
    for chrom_id, chrom_size in chrom_id_size.items():
        seq[chrom_id] = initialise_seq(chrom_size, default_seq_character)

    filter_counter = Counter()
    used_counter = Counter()

    with open(vcf, "r", newline="") as f:
        parsed_lines = parse_lines(f)
        for line in parsed_lines:
            filter_val = line.get("FILTER", "PASS")
            filter_counter[filter_val] += 1

            seq_id = line["CHROM"]
            ref = line["REF"]  # reference base
            alt = line["ALT"]  # alternative base
            pos = parse_position(line["POS"])
            if seq_id in seq and filter_val in ("PASS", "."): # support unfiltered and hard/soft filter
                used_counter[filter_val] += 1
                seq[seq_id][pos - 1] = get_nt_to_add(
                    ref, alt, default_seq_character
                )

    log_filter_info(filter_counter, used_counter)

    for chrom_id, chrom_seq in seq.items():
        assert len(chrom_seq) == chrom_id_size[chrom_id]
        seq[chrom_id] = "".join(chrom_seq)

    return seq


def write_seq(seq: dict, seq_id: str, output_fasta: str = "out.fa") -> None:
    """
    Writes sequences to a FASTA file.

    Parameters:
    seq (dict): A dictionary containing sequence data, where keys are chromosome IDs
        and values are the corresponding sequences.
    seq_id (str): Identifier for the sequence data.
    output_fasta (str, optional): Path to the output FASTA file. Defaults to "out.fa".

    Returns:
    None

    This function writes the sequences from the input dictionary 'seq' to a FASTA file
    specified by 'output_fasta'. Each sequence is written with its corresponding
    chromosome ID as the header, preceded by the 'seq_id'. Sequences are written in the format:

    >{seq_id}_{chrom_id}
    {sequence}

    Example:
    If seq is {'chr1': 'ATCG', 'chr2': 'GCTA'}, seq_id is 'sample', and output_fasta
    is 'output.fasta', the function will write the following content to 'output.fasta':
    >sample_chr1
    ATCG
    >sample_chr2
    GCTA
    """
    with open(output_fasta, "w") as f:
        for chrom_id, sequence in seq.items():
            seq_header = f">{seq_id}_{chrom_id}"
            f.write(f"{seq_header}\n{sequence}\n")


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
        "--default_seq_char",
        "-sc",
        default="N",
        help="Default character for ambiguous bases",
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

    seq = get_seq(args.vcf, args.ref_index, args.default_seq_char)
    write_seq(seq, args.sample, args.output)
