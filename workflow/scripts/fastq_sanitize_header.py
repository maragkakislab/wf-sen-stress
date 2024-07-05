#!/usr/bin/env python

"""
Reads a FASTQ file and sanitizes the header by removing excess information.
"""

import argparse
import sys

def get_file_object(filename):
    """
    Returns a file objects that reads from filename; reads from STDIN if
    filename is -
    """
    if filename == "-":
        return sys.stdin
    return open(filename, "r")

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-f", "--fastq", required=True, help="Input FASTQ file")
args = parser.parse_args()

f = get_file_object(args.fastq)
for ln in f:
    if ln[0] != '@':
        continue

    header = ln.strip().split()
    seq = next(f).strip()
    inter = next(f).strip()
    qual = next(f).strip()

    print(header[0])
    print(seq)
    print(inter)
    print(qual)
