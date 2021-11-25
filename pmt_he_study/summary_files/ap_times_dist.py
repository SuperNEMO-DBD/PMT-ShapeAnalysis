import sys

sys.path.insert(1, '../')

# import ROOT and bash commands
import ROOT

# import python plotting and numpy modules
from format_plot import *
import numpy as np


def main():
    filenames = ["/Users/williamquinn/Desktop/set_4/S95_A25/191106_A1400_B1400_t1613_output.root",
                 "/Users/williamquinn/Desktop/set_4/S95_A25/200213_A1400_B1400_t1037_output.root",
                 "/Users/williamquinn/Desktop/set_4/S95_A25/201106_A1400_B1400_t1200_output.root"]
    exposures = [0, 0.98, 30.28]
    fig = plt.figure(figsize=figsize, facecolor='white')

    markers = ["s", "^", "o"]
    lower, upper = 800, 7000
    n_bins = int((upper - lower)/50)

    for index, i_file in enumerate(filenames):
        file = ROOT.TFile(i_file, "READ")
        date = i_file.split("/")[-1].split("_")[0]
        hist = file.Get(date + "_GAO607_apulse_times_1400V")
        num_hist = file.Get(date + "_GAO607_apulse_num_1400V")
        num = num_hist.GetEntries()
        times = []
        for i in range(hist.GetNbinsX()):
            for j in range(int(hist.GetBinContent(i + 1))):
                times.append(int(hist.GetBinLowEdge(i + 1)))
        file.Close()
        del hist

        freq, bin_edges = np.histogram(times, n_bins, range=(lower, upper))
        width = bin_edges[-1] - bin_edges[-2]
        bin_centres = bin_edges[:-1] + width / 2

        plt.bar(bin_centres, freq/num, width=width, alpha=0.3, color='C{}'.format(index))
        plt.plot(bin_centres, freq/num, "C{}{}".format(index, markers[index]), markersize=marker_size,
                 label='{} atm-day'.format(exposures[index]))

    '''handles, labels = plt.gca().get_legend_handles_labels()
    patch = patches.Patch(color='white', label='Helium Exposure /atm-day')
    handles.insert(0, patch)'''

    plt.xlim(lower, upper)
    plt.yscale('log')
    plt.xlabel('After-Pulse Waveform Time /ns')
    plt.ylabel('Normalised Counts')
    plt.title('After-Pulse Time Distributions')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/ap_times_dist.pdf")


if __name__ == "__main__":
    main()
