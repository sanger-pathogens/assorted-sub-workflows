#!/usr/bin/env python3

## this script is essentially a modified version of the kneaddata functions, check:
# https://github.com/biobakery/kneaddata/blob/fdc9a6bddd8e97446ed6e4d809e680a78bdfc45c/kneaddata/utilities.py#L933 
# https://github.com/biobakery/kneaddata/blob/fdc9a6bddd8e97446ed6e4d809e680a78bdfc45c/kneaddata/run.py#L447
##

import logging
import argparse
import sys

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

def read_file_n_lines(file, n):
    """
    Generator to read a file `n` lines at a time.

    Args:
        file (str): Path to the input file.
        n (int): Number of lines to read at a time. Must be a positive integer.
    
    Yields:
        list[str]: A list of up to `n` lines from the file.

    Raises:
        ValueError: If `n` is not a positive integer.
        FileNotFoundError: If the specified file does not exist.
        IOError: If the file cannot be read.
    """
    if not isinstance(n, int) or n <= 0:
        raise ValueError("`n` must be a positive integer.")

    try:
        with open(file, "r") as file_handle:
            line_batch = []
            for line in file_handle:
                line_batch.append(line)
                if len(line_batch) == n:
                    yield line_batch
                    line_batch = []
            if line_batch:
                yield line_batch
    except FileNotFoundError:
        raise FileNotFoundError(f"Error: File not found: {file}")
    except IOError as e:
        raise IOError(f"Error reading file {file}: {e}")


def remove_repeats_from_fastq(input_fastq, trf_output, output_fastq):
    """
    Remove sequences identified by TRF as containing repeats from a FASTQ file.

    Args:
        input_fastq (str): Path to the input FASTQ file.
        trf_output (str): Path to the TRF output file containing sequence IDs of repeats.
        output_fastq (str): Path to the output FASTQ file (filtered).

    Returns:
        int: The number of sequences removed.
    """
    sequences_with_repeats = set()

    # Read TRF output to collect sequence IDs with repeats
    try:
        c = 0
        with open(trf_output, "r") as file_handle:
            for line in file_handle:
                if line.startswith("@"):  # Sequence IDs in TRF output start with "@"
                    sequences_with_repeats.add(line.strip())
                    c+=1
        logger.info(f" {c} total sequences with repeats from {trf_output}.")
    except FileNotFoundError:
        raise FileNotFoundError(f"Error: TRF output file not found: {trf_output}")
    except Exception as e:
        raise RuntimeError(f"Error reading TRF output file: {e}")

    # Process the FASTQ file and filter sequences
    removed_sequences = 0
    try:
        with open(output_fastq, "w") as file_handle_write:
            for lines in read_file_n_lines(input_fastq, 4):
                sequence_id = lines[0].strip().replace(">", "@")
                # skip reads if present in the TRF file
                if sequence_id in sequences_with_repeats:
                    removed_sequences += 1
                # write read if not
                else:
                    file_handle_write.write("".join(lines))
    except FileNotFoundError:
        raise FileNotFoundError(f"Error: Input FASTQ file not found: {input_fastq}")
    except Exception as e:
        raise RuntimeError(f"Error processing FASTQ file: {e}")

    # Log the result
    logger.info(f"Filtered {removed_sequences} sequences with repeats from {input_fastq} and saved to {output_fastq}.")
    return removed_sequences


def main():
    """
    Main function to parse arguments and call the remove_repeats_from_fastq function.
    """
    parser = argparse.ArgumentParser(
        description="Remove sequences from a FASTQ file that are identified as repeats by TRF."
    )

    # Define arguments
    parser.add_argument(
        "-i", "--input_fastq",
        required=True,
        help="Path to the input FASTQ file."
    )
    parser.add_argument(
        "-t", "--trf_output",
        required=True,
        help="Path to the TRF output file containing sequence IDs with repeats."
    )
    parser.add_argument(
        "-o", "--output_fastq",
        required=True,
        help="Path to the output FASTQ file with filtered sequences."
    )

    args = parser.parse_args()

    # Validate arguments and execute the function
    try:
        removed_count = remove_repeats_from_fastq(args.input_fastq, args.trf_output, args.output_fastq)
        logger.info(f"Successfully removed {removed_count} sequences with repeats.")
    except Exception as e:
        logger.error(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
