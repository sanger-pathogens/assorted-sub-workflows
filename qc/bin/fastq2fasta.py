#!/usr/bin/env python3
import argparse

# This script is essentially a modified version of the function provided in kneaddata, check the link bellow:
# https://github.com/biobakery/kneaddata/blob/master/kneaddata/db_preprocessing/fastq_to_fasta.py#L4

def fastq_to_fasta(fastq_in, fasta_out=None):
    """
    Converts a FASTQ file to a FASTA file.
    
    Args:
        fastq_in (str): Input FASTQ file path.
        fasta_out (str, optional): Output FASTA file path. If not provided, prints to stdout.
    
    Returns:
        bool: True if conversion succeeds, False otherwise.
    """
    try:
        with open(fastq_in, "r") as fp_fastq_in:
            # Open output file if specified
            fp_fasta_out = open(fasta_out, "w") if fasta_out else None
            
            # Process the FASTQ file
            while True:
                # Read four lines at a time (FASTQ record)
                sequence_id = fp_fastq_in.readline().strip()
                if not sequence_id:  # End of file
                    break
                
                sequence = fp_fastq_in.readline().strip()
                plus_line = fp_fastq_in.readline().strip()
                quality_scores = fp_fastq_in.readline().strip()
                
                # Validate FASTQ format (basic check)
                if not sequence_id.startswith("@") or not plus_line.startswith("+"):
                    raise ValueError("Input file does not appear to be in valid FASTQ format.")
                
                # Convert sequence ID and write to output
                fasta_id = sequence_id.replace("@", ">", 1)
                if fp_fasta_out:
                    fp_fasta_out.write(f"{fasta_id}\n{sequence}\n")
                else:
                    print(fasta_id)
                    print(sequence)
            
            # Close output file if opened
            if fp_fasta_out:
                fp_fasta_out.close()
        
        return True  # Success
    
    except FileNotFoundError:
        print(f"Error: File {fastq_in} not found.")
        return False
    except ValueError as e:
        print(f"Error: {e}")
        return False
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return False


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert a FASTQ file to FASTA format.")
    parser.add_argument(
        "fastq_in",
        type=str,
        help="Path to the input FASTQ file."
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
        default=None,
        help="Path to the output FASTA file. If not provided, outputs to stdout."
    )
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Run the conversion
    success = fastq_to_fasta(args.fastq_in, args.output)
    if success:
        print("Conversion completed successfully.")
    else:
        print("Conversion failed.")
        exit(1)
