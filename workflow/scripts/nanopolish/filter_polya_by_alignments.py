#!/usr/bin/env python

"""
Converts the 5' and 3' coordinates of reads into coordinates on meta-features.
In essence, each feature is divided into bins; the 5' and 3' coordinates of
each read contained in each feature are replaced by the bin index in which
they are included. For example for a feature of length 100 containing a read
with coordinates [5', 3']: [15, 95] the corresponding meta-coordinates for
10 bins are [1, 9]. THe script works only for +ve strands and BED coordinates
grater than certain length of nucleotide set in 'min-len' parameter.
"""

import sys
import pysam
import argparse


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-i", "--polya",
                    help="Input nanopolish polya file")
parser.add_argument("-b", "--bam",
                    help="Input SAM/BAM file")
args = parser.parse_args()


aligned_read_ids = {}
infile = pysam.AlignmentFile(args.bam, 'rb')
for read in infile:
    if read.is_unmapped: #check for unmapped and -ve strand
        continue
    qname = read.query_name
    aligned_read_ids[qname] = True


with open(args.polya, 'r') as f:
    print(f.readline(), end='')
    for l in f:
        qname = l.split("\t")[1] 
        if qname in aligned_read_ids:
            print(l, end='')
