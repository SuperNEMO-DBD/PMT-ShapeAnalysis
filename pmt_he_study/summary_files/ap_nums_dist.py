import sys

sys.path.insert(1, '../..')

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

    markers = ["s", "^", "o"]
    lower, upper = 0, 20
    n_bins = 20

    t_nums = []
    t_he_nums = []

    for index, i_file in enumerate(filenames):
        file = ROOT.TFile(i_file, "READ")
        date = i_file.split("/")[-1].split("_")[0]
        num_hist = file.Get(date + "_GAO607_apulse_num_1400V")
        num = num_hist.GetEntries()
        nums = []
        for i in range(num_hist.GetNbinsX()):
            for j in range(int(num_hist.GetBinContent(i + 1))):
                nums.append(int(num_hist.GetBinLowEdge(i + 1)))
        he_num_hist = file.Get(date + "_GAO607_he_apulse_num_1400V")
        he_num = he_num_hist.GetEntries()
        he_nums = []
        for i in range(he_num_hist.GetNbinsX()):
            for j in range(int(he_num_hist.GetBinContent(i + 1))):
                he_nums.append(int(he_num_hist.GetBinLowEdge(i + 1)))

        del he_num_hist
        del num_hist
        file.Close()

        t_nums.append([nums, num])
        t_he_nums.append([he_nums, he_num])

    fig = plt.figure(figsize=figsize, facecolor='white')
    fig.supylabel('Counts per Waveform')
    frame1 = fig.add_axes((.15, .54, .6, .37))
    frame1.set_xticklabels([])

    for i in range(len(t_nums)):
        freq, bin_edges = np.histogram(t_nums[i][0], n_bins, range=(lower, upper))
        width = bin_edges[-1] - bin_edges[-2]
        bin_centres = bin_edges[:-1] + width / 2

        plt.bar(bin_centres, freq/t_nums[i][1], width=width, alpha=0.5, color='C{}'.format(i))
        plt.errorbar(bin_centres, freq/t_nums[i][1], yerr=np.sqrt(freq)/t_nums[i][1], fmt="C{}.".format(i),
                     markersize=1, capsize=cap_size, linewidth=line_width,
                     capthick=cap_thick, label='{:.2f}'.format(exposures[i]))
    plt.yscale('log')
    plt.xlim(lower, upper)

    frame2 = fig.add_axes((.15, .14, .6, 0.37))
    for i in range(len(t_he_nums)):
        freq, bin_edges = np.histogram(t_he_nums[i][0], n_bins, range=(lower, upper))
        width = bin_edges[-1] - bin_edges[-2]
        bin_centres = bin_edges[:-1] + width / 2

        plt.bar(bin_centres, freq / t_he_nums[i][1], width=width, alpha=0.5, color='C{}'.format(i))
        plt.errorbar(bin_centres, freq / t_he_nums[i][1], yerr=np.sqrt(freq) / t_he_nums[i][1], fmt="C{}.".format(i),
                     markersize=1, capsize=cap_size, linewidth=line_width,
                     capthick=cap_thick, label='{:.2f}'.format(exposures[i]))

    plt.xlim(lower, upper)
    plt.yscale('log')
    plt.xlabel('After-Pulse Number')

    handles, labels = plt.gca().get_legend_handles_labels()
    patch = patches.Patch(color='white', label='Exposure')
    patch1 = patches.Patch(color='white', label='atm-day')
    handles.insert(0, patch1)
    handles.insert(0, patch)
    plt.legend(handles=handles, loc='best', bbox_to_anchor=(1.4, 2))

    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/ap_nums_dist.pdf")


if __name__ == "__main__":
    main()
