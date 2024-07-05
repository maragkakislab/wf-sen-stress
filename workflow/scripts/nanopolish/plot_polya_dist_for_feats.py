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


def parse_args(args):

    # Create an ArgumentParser object
    parser = argparse.ArgumentParser()

    # Add a positional argument for the filename
    parser.add_argument('ifile',
                        help='Name of the nanopolish polya tab-delimited file')
    parser.add_argument('ofile',
                        help='Name of output file with figure')
    parser.add_argument("-b","--bed",
                        help = "BED file with features.")
    parser.add_argument("-m","--bam",
                        help = "BAM file with alignments. Must be indexed.")
    return parser.parse_args()


def main():
    # Parse the command-line arguments
    args = parse_args(sys.argv[1:])

    # Read the file with poly(A) tails into a pandas DataFrame
    df = pd.read_csv(args.ifile, delimiter='\t')

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


        num_of_plots = 1
        inch_per_plot = 4 # for pdf it essentially only helps to scale text
        fig, axes = plt.subplots(1, num_of_plots,
                                 figsize=(inch_per_plot*num_of_plots, inch_per_plot),
                                 layout='constrained')
        sns.ecdfplot(ax=axes, data=df, x="polya_length", hue="feat")
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
    # Get the values of the histogram
    hrange = (0, 256)
    values, bin_edges = np.histogram(df['polya_length'], bins=hrange[1]-hrange[0], range=hrange)
    values = scipy.stats.zscore(values)

    # Perform the Fourier transformation
    data_fft = np.fft.fft(values)

    # Calculate the absolute value of the complex numbers
    data_fft_abs = np.abs(data_fft)

    # Calculate the Power Spectral Density
    data_fft_psd = np.abs(data_fft)**2
    data_fft_db = 10*np.log10(data_fft_psd)

    # Get the frequency values
    frequency = np.fft.fftfreq(len(values))
    positives = frequency > 0 # For real signal negative frequencies are not required

    # Calculate the period
    period = 1 / frequency

    # Use the Welch method
    wf, wpsd = welch(values, fs=1, nperseg=32, noverlap=16)

    # Setup the figure.
    num_of_plots = 7
    inch_per_plot = 4 # for pdf it essentially only helps to scale text
    fig, axes = plt.subplots(1, num_of_plots,
                             figsize=(inch_per_plot*num_of_plots, inch_per_plot),
                             layout='constrained')
    axes_it = iter(axes) # An iterator for convenient looping on the axes.

    # Create plots
    ax = next(axes_it)
    summary= df[['polya_length']].describe()
    ax.text(0.1, 0.5, str(summary))
    ax.set_title('Poly(A) tail length statistics')

    # Plot the poly(A) tail length distribution
    ax = next(axes_it)
    sns.histplot(ax=ax, data=df, x='polya_length')
    ax.set_xlim(0,500)
    ax.set_title('Poly(A) tail length distribution')

    # Plot the spectra (frequency)
    ax = next(axes_it)
    ax.plot(frequency[positives], data_fft_db[positives])
    ax.set_title('Fourier Transform Spectra')
    ax.set_xlabel('Frequency')
    ax.set_ylabel('Amplitude (db)')

    # Plot the spectra (period)
    ax = next(axes_it)
    ax.plot(period[positives], data_fft_db[positives])
    ax.set_title('Fourier Transform Spectra')
    ax.set_xlabel('Period (nts)')
    ax.set_ylabel('Amplitude (db)')

    # Zoom-in (period)
    ax = next(axes_it)
    ax.plot(period[positives], data_fft_db[positives])
    ax.set_xlim(0,20)
    ax.set_title('Fourier Transform Spectra (zoomed-in)')
    ax.set_xlabel('Period (nts)')
    ax.set_ylabel('Amplitude (db)')

    # Use the Welch method.
    ax = next(axes_it)
    ax.plot(wf, wpsd)
    ax.set_title('Welch Transform')
    ax.set_xlabel('Frequency')
    ax.set_ylabel('Power Spectral Density')

    # Use the Welch method (zoom-in the period).
    ax = next(axes_it)
    ax.plot(1/wf, wpsd)
    ax.set_xlim(0,20)
    ax.set_title('Welch Transform (zoomed-in)')
    ax.set_xlabel('Period (nts)')
    ax.set_ylabel('Power Spectral Density')

    return fig


if __name__ == "__main__":
    main()
