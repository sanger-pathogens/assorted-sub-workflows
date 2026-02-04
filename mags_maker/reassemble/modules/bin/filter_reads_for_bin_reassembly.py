#!/usr/bin/env python3

import os
import argparse
import gzip

from pysam import AlignmentFile

complement = {
    'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
    'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
    'N': 'N', 'n': 'n'
}

def rev_comp(seq):
    rev_comp = ""
    for n in seq:
        rev_comp += complement[n]
    return rev_comp[::-1]

def main():
    parser = argparse.ArgumentParser(
        description="Filter mapped reads for bin reassembly based on mismatch thresholds."
    )
    parser.add_argument("original_bin_folder", help="Folder containing original bin fasta files")
    parser.add_argument("output_dir", help="Directory to store output FASTQ files")
    parser.add_argument("strict_snp_cutoff", type=int, help="Strict mismatch threshold")
    parser.add_argument("permissive_snp_cutoff", type=int, help="Permissive mismatch threshold")
    parser.add_argument(
        "--mapped_reads", type=str, default=None,
        help="Mapped reads (sam/bam) input instead of stdin. If omitted, reads are expected from stdin."
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    strict_snp_cutoff = args.strict_snp_cutoff
    permissive_snp_cutoff = args.permissive_snp_cutoff

    print("loading contig to bin mappings...")
    contig_bins = {}
    for bin_file in os.listdir(args.original_bin_folder):
        if bin_file.endswith(".fa") or bin_file.endswith(".fasta"):
            bin_name = ".".join(bin_file.split("/")[-1].split(".")[:-1])
            with open(os.path.join(args.original_bin_folder, bin_file)) as f:
                for line in f:
                    if line[0] != ">":
                        continue
                    contig_bins[line[1:-1]] = bin_name

    print("Parsing sam file and writing reads to appropriate files depending what bin they alligned to...")
    files = {}
    opened_bins = {}

    if args.mapped_reads.endswith(".bam"):
        input_stream = AlignmentFile(args.mapped_reads, "rb")
    elif args.mapped_reads.endswith(".sam"):
        input_stream = AlignmentFile(args.mapped_reads, "r")
    else:
        input_stream = AlignmentFile("-", "r")
    
    # Includes unmapped reads (probably don't want this)
    # 
    sam_iter = input_stream.fetch(until_eof=True)

    for read in sam_iter:
        line = read.to_string()
        cut = line.strip().split("\t")
        binary_flag = bin(int(cut[1]))

        if binary_flag[-7] == "1":
            F_line = line
            continue
        elif binary_flag[-8] == "1":
            R_line = line

            F_cut = F_line.strip().split("\t")
            R_cut = R_line.strip().split("\t")

            if F_cut[2] == "*" and R_cut[2] == "*":
                continue

            if F_cut[2] != R_cut[2]:
                if F_cut[2] not in contig_bins or R_cut[2] not in contig_bins:
                    continue
                bin1 = contig_bins[F_cut[2]]
                bin2 = contig_bins[R_cut[2]]
                if bin1 != bin2:
                    continue
                bin_name = bin1
            else:
                contig = F_cut[2]
                if contig not in contig_bins:
                    continue
                bin_name = contig_bins[contig]

            if "NM:i:" not in F_line and "NM:i:" not in R_line:
                continue

            if bin_name not in opened_bins:
                opened_bins[bin_name] = None
                for mode in ["strict", "permissive"]:
                    for end in [1, 2]:
                        filename = os.path.join(args.output_dir, f"{bin_name}.{mode}_{end}.fastq.gz")
                        files[filename] = gzip.open(filename, "wt")

            cumulative_mismatches = 0
            for field in F_cut:
                if field.startswith("NM:i:"):
                    cumulative_mismatches += int(field.split(":")[-1])
                    break
            for field in R_cut:
                if field.startswith("NM:i:"):
                    cumulative_mismatches += int(field.split(":")[-1])
                    break

            F_binary_flag = bin(int(F_cut[1]))
            R_binary_flag = bin(int(R_cut[1]))

            if F_binary_flag[-5] == '1':
                F_cut[9] = rev_comp(F_cut[9])
                F_cut[10] = F_cut[10][::-1]
            if R_binary_flag[-5] == '1':
                R_cut[9] = rev_comp(R_cut[9])
                R_cut[10] = R_cut[10][::-1]

            if cumulative_mismatches < strict_snp_cutoff:
                files[os.path.join(args.output_dir, f"{bin_name}.strict_1.fastq.gz")].write(f"@{F_cut[0]}/1\n{F_cut[9]}\n+\n{F_cut[10]}\n")
                files[os.path.join(args.output_dir, f"{bin_name}.strict_2.fastq.gz")].write(f"@{R_cut[0]}/2\n{R_cut[9]}\n+\n{R_cut[10]}\n")

            if cumulative_mismatches < permissive_snp_cutoff:
                files[os.path.join(args.output_dir, f"{bin_name}.permissive_1.fastq.gz")].write(f"@{F_cut[0]}/1\n{F_cut[9]}\n+\n{F_cut[10]}\n")
                files[os.path.join(args.output_dir, f"{bin_name}.permissive_2.fastq.gz")].write(f"@{R_cut[0]}/2\n{R_cut[9]}\n+\n{R_cut[10]}\n")

    print("closing files")
    if args.sam:
        input_stream.close()
        
    for f in files.values():
        f.close()

    print("Finished splitting reads!")

if __name__ == "__main__":
    main()
