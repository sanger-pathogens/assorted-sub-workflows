#!/usr/bin/env python3
import sys


class ContigParsingWarning(Exception):
    def __init__(self, line):
        msg = f"Warning: Couldn't parse integer from line: {line.strip()}\n"
        print(msg, file=sys.stderr)


class NoValidContigsWarning(Exception):
    def __init__(self, filename, min_len):
        msg = f"Warning: No contigs found longer than {min_len} in {filename}\n"
        print(msg, file=sys.stderr)


class InvalidArgumentsError(Exception):
    def __init__(self, expected, given):
        msg = f"Expected {expected} arguments, but got {given}.\n"
        print(msg, file=sys.stderr)


try:
    if len(sys.argv) != 3:
        raise InvalidArgumentsError(2, len(sys.argv) - 1)
except IOError as ex:
    print(f"{ex.args[0]}", file=sys.stderr)
    sys.exit(1)
except Exception as ex:
    print(f"Unexpected error: {ex.args[0]}")


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
                ContigParsingWarning(line)
                continue

            if value < min_len:
                break
            else:
                print(line.strip())
                print("working")
                contigs_found = True

if not contigs_found:
    NoValidContigsWarning(filename, min_len)
    # Indicate that no valid contigs were found
