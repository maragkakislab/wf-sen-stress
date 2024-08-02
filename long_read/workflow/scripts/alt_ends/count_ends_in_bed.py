"""
Plots the distribution of reads in a BAM file around a given set of features
in a BED file. The user can specify which end (5' or 3') of the reads and the
features will be used as reference for the comparison. For example: We assume
that the user selects the 5' end of the reads and the 5' end of the features
as reference. Then a read that maps at position 10 of chr1 will be at a
relative position of -5 nt compared to a feature aligning at position 15 of
chr1. The same concept is applied for all reads against all features and a
distribution of relative positions is constructed.
"""


import pysam
import argparse


class Feature:
    def __init__(self, chrom, start, end, name, score, strand):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.name = name
        self.score = score
        self.strand = strand

    def __str__(self):
        return "{}\t{}\t{}\t{}\t{}\t{}".format(
                                               self.chrom, self.start,
                                               self.end, self.name,
                                               self.score, self.strand)


def find_pos(start, end, strand, pos):
    if pos != '5p' and pos != '3p':
        raise ValueError('Incorrectly specified position')
    if strand != '+' and strand != '-':
        raise ValueError('Incorrectly specified strand')
    final_pos = start
    if strand == '-' and pos == '5p' or strand == '+' and pos == '3p':
        final_pos = end
    return final_pos


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-m","--bam", required=True,
                    help = "BAM file with reads. Must be indexed.")
    parser.add_argument("-b","--bed", required=True,
                    help = "BED file with features.")
    parser.add_argument("--extend-start", type = int, default = 0,
                    help = "Number of nucleotides to extend bed start (default: %(default)s).")
    parser.add_argument("--extend-end", type = int, default = 0,
                    help = "Number of nucleotides to extend bed end (default: %(default)s).")
    parser.add_argument("-r","--rpos", default='5p',
                    help = "Reference point for reads; one of 5p or 3p (default: %(default)s)")
    args = parser.parse_args()

    # Open the bam file
    bamfile = pysam.AlignmentFile(args.bam, "rb")

    # Loop on the BED file and query the BAM file to get overlapping reads.
    feature_count = 0
    with open(args.bed) as bed:
        for line in bed:
            feature_count += 1
            cols = line.strip().split('\t')
            feat = Feature(cols[0], int(cols[1]), int(cols[2]), cols[3],
                           cols[4], cols[5])
            feat.score = 0 # reset score to use for counting overlapping reads

            if bamfile.header.get_tid(feat.chrom) == -1:
                continue

            # define region of interest
            istart = feat.start - args.extend_start
            if istart < 0:
                istart = 0
            iend = feat.end + args.extend_end
            # TODO check if iend is outside the chromosome bounds

            reads = bamfile.fetch(feat.chrom, istart, iend)
            for r in reads:
                rs = '+' if r.is_forward else '-'
                pos = find_pos(r.reference_start, r.reference_end - 1, rs,
                               args.rpos)
                if pos >= istart and pos <= iend - 1:
                    feat.score += 1
            print(feat)


if __name__ == "__main__":
    main()
