import argparse
import scipy
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.signal import welch

# Create an ArgumentParser object
parser = argparse.ArgumentParser()

# Add a positional argument for the filename
parser.add_argument('ifile', help='Name of the nanopolish polya tab-delimited file')
parser.add_argument('ofile', help='Name of output file with figure')

# Parse the command-line arguments
args = parser.parse_args()

# Read the file into a pandas DataFrame
df = pd.read_csv(args.ifile, delimiter='\t')

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

f, psd = welch(values, fs=1, nperseg=32, noverlap=16)

# Use Seaborn to plot a histogram of the length column
with PdfPages(args.ofile) as pages:
    sns.histplot(data=df, x='polya_length')
    pages.savefig()
    plt.close()

    sns.kdeplot(df['polya_length'])
    pages.savefig()
    plt.close()

    # Plot the spectra (frequency)
    plt.plot(frequency[positives], data_fft_db[positives])
    plt.title('Fourier Transform Spectra')
    plt.xlabel('Frequency')
    plt.ylabel('Amplitude')
    pages.savefig()
    plt.close()

    # Plot the spectra (period)
    plt.plot(period[positives], data_fft_db[positives])
    plt.title('Fourier Transform Spectra')
    plt.xlabel('Period (nts)')
    plt.ylabel('Amplitude (db)')
    pages.savefig()
    plt.close()

    # Zoom-in (period)
    plt.plot(period[positives], data_fft_db[positives])
    plt.xlim(0,20)
    plt.title('Fourier Transform Spectra (zoom-in)')
    plt.xlabel('Period (nts)')
    plt.ylabel('Amplitude (db)')
    pages.savefig()
    plt.close()


    plt.plot(1/f, psd, 'o-')
    plt.xlabel('Period (nts)')
    plt.ylabel('Power Spectral Density')
    pages.savefig()
    plt.close()

    plt.plot(1/f, psd)
    plt.xlim(0,20)
    plt.ylim(0,0.1)
    plt.xlabel('Period (nts)')
    plt.ylabel('Power Spectral Density')
    pages.savefig()
    plt.close()
