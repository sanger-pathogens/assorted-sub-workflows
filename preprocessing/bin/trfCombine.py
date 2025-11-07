#!/usr/bin/env python3
import argparse

def read_trf(file_path):
    """Read a TRF file and return a dictionary of read identifiers with their corresponding lines."""
    reads = {}
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('@'):  # Read identifier lines start with '@'
                read_id = line.strip()
                # Remove the trailing number (1 or 2) from the read identifier
                read_id_base = read_id.rsplit('/', 1)[0]  # Get the part before the slash
                reads[read_id_base] = read_id  # Store the original read ID
    return reads

def process_read(reads_source, reads_target, read_id_base, output_file, unpaired_reads_fl, suffix_replace):
    """
    Process a read from one source, attempt to find its equivalent in the target, 
    and handle unpaired reads.
    
    Args:
        reads_source (dict): The source reads dictionary (e.g., reads_trf1 or reads_trf2).
        reads_target (dict): The target reads dictionary to check for equivalence.
        read_id_base (str): The base read ID being processed.
        output_file (file object): File object for writing combined reads.
        unpaired_reads_fl (file object): File object for writing unpaired reads.
        unpaired_reads_lst (list): List to collect unpaired reads.
        suffix_replace (tuple): Tuple of (source_suffix, target_suffix) for ID replacement.
    """
    try:
        # Write the primary read to the output file
        output_file.write(reads_source[read_id_base] + '\n')

        # Generate the equivalent read ID by replacing suffixes
        source_suffix, target_suffix = suffix_replace

        # sanity check read format
        assert(len(reads_source[read_id_base].split("/")) == 2)

        equivalent_suffix = reads_source[read_id_base].split("/")[1].replace(source_suffix, target_suffix)
        equivalent_read = f"{reads_source[read_id_base].split('/')[0]}/{equivalent_suffix}"

        # Check if the equivalent read exists in the target
        if equivalent_read not in reads_target:
            # store on the combined trf (reads to be removed)
            output_file.write(equivalent_read + '\n')
            # store on the unpaired_reads (for the record)
            unpaired_reads_fl.write(equivalent_read + '\n')

    except AssertionError:
        print(f"ERROR: {reads_source[read_id_base]} is not a valid read ID.")
        print(f"       example of a valid read format: '@A01404:579:HVVNFDRX5:2:2122:28076:11569/1'")
        exit(1)

def generate_combined_trf(trf1_path, trf2_path, output_path, unpaired_path):
    """Generate a new TRF file that includes equivalent reads from both TRF files."""
    reads_trf1 = read_trf(trf1_path)
    reads_trf2 = read_trf(trf2_path)
    # store reads which lose their pairs
    unpaired_reads_fl = open(unpaired_path,"w")
    with open(output_path, 'w') as output_file:
        #  get list of unique identifiers 
        trfs_reads_set = list(set(reads_trf1.keys()).union(reads_trf2.keys()))
        # sort list so output always have same order
        # NOTE: specially usefull for testing the code
        trfs_reads_set.sort()
        for read_id_base in trfs_reads_set:
            if read_id_base in reads_trf1:
                process_read(
                    reads_source=reads_trf1,
                    reads_target=reads_trf2,
                    read_id_base=read_id_base,
                    output_file=output_file,
                    unpaired_reads_fl=unpaired_reads_fl,
                    suffix_replace=("1", "2")
                )
                continue

            if read_id_base in reads_trf2:
                process_read(
                    reads_source=reads_trf2,
                    reads_target=reads_trf1,
                    read_id_base=read_id_base,
                    output_file=output_file,
                    unpaired_reads_fl=unpaired_reads_fl,
                    suffix_replace=("2", "1")
                )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a combined TRF file from two TRF files.")
    parser.add_argument("-f1", "--trf1", required=True, help="Path to the first TRF file.")
    parser.add_argument("-f2", "--trf2", required=True, help="Path to the second TRF file.")
    parser.add_argument("-o", "--output", required=True, help="Output path for the combined TRF file.")
    parser.add_argument("-u", "--unpaired", required=True, help="Output path for the unpaired reads ids list.")

    args = parser.parse_args()

    generate_combined_trf(args.trf1, args.trf2, args.output, args.unpaired)

    print(f"Combined TRF file created at: {args.output}")