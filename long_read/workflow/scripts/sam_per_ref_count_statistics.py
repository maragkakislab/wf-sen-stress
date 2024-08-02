#!/usr/bin/env python

import argparse
import statistics as st
import pysam

"""
Script to calculate the reference counts and length statistics (like read
count, mean length, median length, stdev, variance, min length, max length,
 total reads and mapped reads) of BAM/SAM file.
"""

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--ifile",
    help="Seq file name with read information as input in sam/bam format")
parser.add_argument("-f", "--sam", action='store_true',
    help="Input file format only if SAM file; default BAM format")
parser.add_argument("-r", "--ref-col-name",
    default="reference", help="Reference output column name, default:reference")
parser.add_argument("-c", "--cnt-col-name",
    default="count", help="Read count output column name, default: count")
parser.add_argument("-n", "--opt-col-name",
    help="Name of an optional column e.g. sample_name")
parser.add_argument("-v", "--opt-col-val",
    help="Value for the optional column; same for all rows")
parser.add_argument("-s", "--col_delimiter", default="\t",
    help="Delimiter to seperate the columns of the output file, default : tab(\t)")
parser.add_argument("-x", "--no-zeros",
    action='store_true', help="Skip references with 0 reads")
args = parser.parse_args()

#Arguments about input file type, BAM format is default
ifiletype = "rb"
if args.sam:
    ifiletype = "r"

#Creating the dictionary for reads
reference_counts = {}

#Initiating -- total and mapped reads
total_reads = 0
mapped_reads = 0


#Opening and processing the input file
bamfile = pysam.AlignmentFile(args.ifile, ifiletype)

#Processing data for the reads (mapped)
for seq in bamfile:
    total_reads += 1
    if seq.is_unmapped:
        continue
    mapped_reads += 1
    reference = seq.reference_name
    ref_len = seq.reference_length
    query_length = seq.query_length
    if reference not in reference_counts:
        reference_counts[reference] = []
    reference_counts[reference] += [query_length]

# If transcripts with 0 counts are not explicitly excluded, then add them to
# the dictionary. Assumes that the SAM header exists.
if not args.no_zeros:
    header = bamfile.header
    if 'SQ' in header:
        for elem in header['SQ']:
            ref = elem['SN']
            if ref not in reference_counts:
                reference_counts[ref] = []

delim = args.col_delimiter

#print the header of file
header = [args.ref_col_name,
          args.cnt_col_name,
          "mean_length",
          "median_length",
          "stdev_length",
          "var_length",
          "max_length",
          "min_length",
          # "mode_length"
          "total_reads",
          "mapped_reads"]

if args.opt_col_name and args.opt_col_val:
    header += [args.opt_col_name]
print(delim.join(header))

#print the value per read
for ref, count in reference_counts.items():
    if len(count) == 1:
        row = [ref,
               str(len(count)),
               str(st.mean(count)),
               str(st.median(count)),
               str('inf'),
               str('inf'),
               str(max(count)),
               str(min(count)),
               # str(multimode(count)),
               str(total_reads),
               str(mapped_reads)]
    elif len(count) > 1:
        row = [ref,
               str(len(count)),
               str(st.mean(count)),
               str(st.median(count)),
               str(st.stdev(count)),
               str(st.variance(count)),
               str(max(count)),
               str(min(count)),
               # str(multimode(count)),
               str(total_reads),
               str(mapped_reads)]
    else:
        row = [ref,
               str(len(count)),
               str('NaN'),
               str('NaN'),
               str('NaN'),
               str('NaN'),
               str('NaN'),
               str('NaN'),
               # str('NaN'),
               str(total_reads),
               str(mapped_reads)]

    if args.opt_col_name and args.opt_col_val:
        row += [args.opt_col_val]
    print(delim.join(row))
