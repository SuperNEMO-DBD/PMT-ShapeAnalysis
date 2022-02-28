import sys

sys.path.insert(1, '../../')

# import ROOT and bash commands
import ROOT

# import python plotting and numpy modules
from pmt_he_study.format_plot import *
from functions.other_functions import process_exposure
import numpy as np


def main():
    filenames = ["/Users/williamquinn/Desktop/data/1400V/200306_A1400_B1400_t2134_output.root",
                 "/Users/williamquinn/Desktop/data/1400V/200506_A1400_B1400_t2247_output.root",
                 "/Users/williamquinn/Desktop/data/1400V/200706_A1400_B1400_t0012_output.root",
                 "/Users/williamquinn/Desktop/data/1400V/200906_A1400_B1400_t2242_output.root",
                 "/Users/williamquinn/Desktop/data/1400V/201106_A1400_B1400_t1200_output.root",
                 "/Users/williamquinn/Desktop/data/1400V/211209_A1400_B1400_t1101_output.root"]
    # exposures = process_exposure(np.array([int(filenames[i].split("/")[-1].split("_")[0]) for i in range(len(filenames))]))
    exposures = [3.47, 9.57, 15.67, 21.77, 27.87, 67.58]
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

        #plt.bar(bin_centres, freq/num, width=width, alpha=0.3, color='C{}'.format(index))
        plt.plot(bin_centres, freq/num, "C{}-".format(index), markersize=marker_size,
                 label='{:.2f}'.format(exposures[index]))

    plt.xlim(lower, upper)
    #plt.yscale('log')
    plt.xlabel('After-Pulse Waveform Time /ns')
    plt.ylabel('Counts per Waveform')
    plt.title('After-Pulse Time Distributions')
    handles, labels = plt.gca().get_legend_handles_labels()
    patch = patches.Patch(color='white', label='Exposure')
    patch1 = patches.Patch(color='white', label='atm-day')
    handles.insert(0, patch1)
    handles.insert(0, patch)
    plt.legend(handles=handles, loc='best')
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/ap_times_dist.pdf")


if __name__ == "__main__":
    main()
