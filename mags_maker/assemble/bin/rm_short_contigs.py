#!/usr/bin/env python3
import sys

try:
    if len(sys.argv) != 3:
        raise IOError(f"Usage: {sys.argv[0]} <min_len> <filename>\nNumber of args expected {2} given {len(sys.argv)-1}")
except IOError as ex: 
    print(f"{ex.args[0]}", file=sys.stderr)
    sys.exit(1)
except Exception as ex:
    print(f"Unexpected error: {ex.args[0]}", file=sys.stderr)
    sys.exit(1)


min_len = int(sys.argv[1])
filename = sys.argv[2]

with open(filename) as f:
    contigs_found = False
    for line in f:
        if not line.startswith(">"):
            print(line.strip())
        else:
            try:
                value = int(line.split("_")[3])
            except (IndexError, ValueError):
                print(f"Warning: Couldn't parse integer from line: {line.strip()}")
                continue

            if value < min_len:
                break
            else:
                print(line.strip())
                contigs_found = True

if not contigs_found:
    print(f"Warning: No contigs found longer than {min_len} in {filename}", file=sys.stderr)
    #sys.exit(1)  # Indicate that no valid contigs were found
