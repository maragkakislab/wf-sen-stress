import sys
import argparse
import scipy
import re
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.signal import welch
import pysam


def build_read_index_name_dict(filename):
    my_dict = {}

    with open(filename, "r") as f:
        next(f)  # Skip the header line
        for line in f:
            columns = line.strip().split("\t")
            my_dict[columns[0]] = columns[1]


def parse_args(args):

    # Create an ArgumentParser object
    parser = argparse.ArgumentParser()

    # Add a positional argument for the filename
    parser.add_argument('ifile',
                        help='Name of the nanopolish polya tab-delimited file')
    parser.add_argument('ofile',
                        help='Name of output file with figure')
    parser.add_argument("-s","--summary",
                        help = "Nanopolish summary file with read index to name associations.")
    parser.add_argument("-p","--indiv-prob",
                        help = "m6anet output file with individual probabilities.")
    parser.add_argument("-b","--bed",
                        help = "BED file with features.")
    parser.add_argument("-m","--bam",
                        help = "BAM file with alignments. Must be indexed.")
    return parser.parse_args()


def main():
    # Parse the command-line arguments
    args = parse_args(sys.argv[1:])

    # Read file that associates read index with name.
    df_read_idx_to_name = pd.read_csv(args.summary, sep="\t",
                                      usecols=["read_index", "read_name"])

    # Read file with individual modification probabilities.
    df_indiv_prob = pd.read_csv(args.indiv_prob)

    # Merge the two dataframes.
    df_merged = pd.merge(df_indiv_prob, df_read_idx_to_name,
                         on="read_index")

    # Aggregate probabilities per read_name
    df_merged = df_merged.groupby('read_name')['probability_modified'].max().reset_index()

    # Read the file with poly(A) tails into a pandas DataFrame
    df_polya = pd.read_csv(args.ifile, delimiter='\t')

    # Merge the individual modification probabilities to the polyA lengths.
    df = pd.merge(df_polya, df_merged,
                  left_on="readname", right_on="read_name")

    # Bin modification probabilities into equally sized bins.
    bin_edges = [i/10 for i in range(6)]  # Generate equally spaced values
    bin_edges[-1] = 1  # Set the last bin edge to infinity
    df['mod_group'] = pd.cut(df['probability_modified'], bins=bin_edges,
                             labels=None, include_lowest=True)

    fig_for_all = plot_df(df)
    fig_for_all.suptitle("All")

    figs = []
    if args.bed != None and args.bam != None:
        annotate_df_with_overlapping_feats(df, args.bed, args.bam)
        figs = process_individual_features(df, args.bed)

    with PdfPages(args.ofile) as pages:
        for f in [fig_for_all] + figs:
            pages.savefig(f)
            plt.close()


        if args.bed != None and args.bam != None:
            num_of_plots = 1
            inch_per_plot = 12 # for pdf it essentially only helps to scale text
            fig, axes = plt.subplots(1, num_of_plots,
                                     figsize=(inch_per_plot*num_of_plots, inch_per_plot),
                                     layout='constrained')
            sns.ecdfplot(ax=axes, data=df, x="polya_length",
                         hue=df[['feat', 'mod_group']].apply(tuple, axis=1))
            axes.set_xlim(0,500)
            pages.savefig(fig)
        plt.close()


def annotate_df_with_overlapping_feats(df, bed, bam):
    # Open the bam file
    bamfile = pysam.AlignmentFile(bam, "rb")

    # Loop on BED features.
    with open(bed) as bed_file:
        df['feat'] = None
        for line in bed_file:
            # Region of interest.
            cols = line.strip().split('\t')
            chrom = cols[0]
            start = int(cols[1])
            end = int(cols[2])
            info = cols[9]
            feat_id = re.search(r'ID=([^;]+)', info).group(1)

            # Fetch read ids that are fully contained in region of interest.
            reads = bamfile.fetch(chrom, start, end)
            ids = [r.query_name for r in list(reads) if start <= r.reference_start and r.reference_end <= end]

            # Annotate data frame with new column that indicates the
            # overlapping feat.
            df.loc[df['readname'].isin(ids), 'feat'] = feat_id


def process_individual_features(df, bed):
    figs = [] # will contain one figure for every feature in the BED

    # Loop on BED features.
    with open(bed) as bed_file:
        for line in bed_file:
            # Region of interest.
            cols = line.strip().split('\t')
            chrom = cols[0]
            start = int(cols[1])
            end = int(cols[2])
            info = cols[9]
            feat_id = re.search(r'ID=([^;]+)', info).group(1)

            # Filter data frame to only keep reads in region.
            fdf = df[df['feat'] == feat_id]

            # Make plots for filtered data frame.
            fig = plot_df(fdf)
            fig.suptitle(chrom + ':' + str(start) + '-' + str(end))
            figs.append(fig)

        return figs


def plot_df(df):
    # Setup the figure.
    num_of_plots = 4
    inch_per_plot = 5 # for pdf it essentially only helps to scale text
    fig, axes = plt.subplots(1, num_of_plots,
                             figsize=(inch_per_plot*num_of_plots, inch_per_plot),
                             layout='constrained')
    axes_it = iter(axes) # An iterator for convenient looping on the axes.

    # Create plots

    # Plot the poly(A) tail length distribution
    ax = next(axes_it)
    sns.histplot(ax=ax, data=df, x='polya_length', hue='mod_group')
    ax.set_xlim(0,500)
    ax.set_title('Poly(A) tail length distribution')

    # Plot the poly(A) tail length distribution
    ax = next(axes_it)
    sns.kdeplot(ax=ax, data=df, x='polya_length', hue='mod_group')
    ax.set_xlim(0,500)
    ax.set_title('Poly(A) tail length distribution')

    # Plot the poly(A) tail length cumulative distribution
    ax = next(axes_it)
    sns.ecdfplot(ax=ax, data=df, x="polya_length", hue="mod_group")
    ax.set_xlim(0,500)
    ax.set_title('Poly(A) tail length cum. distribution')

    # Print a summary
    ax = next(axes_it)
    summary= df.groupby('mod_group')['polya_length'].describe()
    ax.text(0.1, 0.5, str(summary))
    ax.set_title('Poly(A) tail length statistics')

    return fig


if __name__ == "__main__":
    main()

