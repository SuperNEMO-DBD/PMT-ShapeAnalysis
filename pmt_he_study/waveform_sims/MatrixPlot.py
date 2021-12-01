import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import ROOT
import sys
sys.path.insert(1, '../')

from format_plot import *
from scipy.optimize import curve_fit

sys.path.insert(1, '../..')
from functions.other_functions import *


def test_model(x, q, m, c):
    return q * x * x + m * x + c


def adjust_pars(M):
    averages = []
    stds = []
    nums = np.array([i for i in range(len(M[0]))])
    for i_row in range(len(M)):
        weights = []
        for i_col in range(len(M[i_row])):
            weights.append(M[i_row][i_col])
        weights = np.array(weights)
        average = np.average(nums, weights=weights)

        std = np.sqrt(np.sum(weights * (nums - average) ** 2) / np.sum(weights))
        stds.append(std)
        averages.append(average)

    means = np.linspace(0, 18, 19*5)
    f = []
    for i_mean in means:
        probs = []
        for n in nums:
            prob = np.exp(-1 * i_mean) * (i_mean ** n) / np.math.factorial(n) * averages[n]
            probs.append(prob)
        f.append(np.sum(probs))

    f = np.array(f)
    # print(len(f))
    upper = 5
    my_list = np.where(means < upper)[0]

    t_popt, t_pcov = curve_fit(f=test_model, xdata=f[my_list], ydata=means[my_list])

    fig1 = plt.figure(figsize=figsize)
    frame1 = fig1.add_axes((.13, .35, .8, .6))
    plt.plot(f[my_list], means[my_list], 'o', markersize=marker_size)
    plt.plot(f[my_list], test_model(f[my_list], *t_popt), label=r'Model $y=qx^2 +mx + c$')
    frame1.set_xticklabels([])
    plt.ylabel("Corrected Average")
    plt.title("Average After-Pulse Number Corrected vs. Measured")
    plt.ylim(means[0], upper)
    plt.xlim(0, 4)
    handles, labels = plt.gca().get_legend_handles_labels()
    patch = patches.Patch(color='white', label=r'$q =$ {:.3f} ± {:.3f}'.format(t_popt[0], np.sqrt(t_pcov[0, 0])))
    patch_1 = patches.Patch(color='white', label=r'$m =$ {:.3f} ± {:.3f}'.format(t_popt[1], np.sqrt(t_pcov[1, 1])))
    patch_2 = patches.Patch(color='white', label=r'$c =$ {:.3f} ± {:.3f}'.format(t_popt[2], np.sqrt(t_pcov[2, 2])))
    handles.extend([patch, patch_1, patch_2])
    plt.legend(handles=handles, loc='best')
    frame2 = fig1.add_axes((.13, .15, .8, .2))
    plt.plot(f[my_list], (means[my_list] - test_model(f[my_list], *t_popt))/test_model(f[my_list], *t_popt),
             'ko', markersize=marker_size)
    # plt.ylim(means[0], upper)
    plt.xlim(0, 4)
    plt.ylim(-0.015, 0.015)
    plt.xlabel("Measured Average")
    plt.ylabel('(data - model)/model', fontsize=5)
    plt.axhline(0, ls='--', color='black')
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/PMT_project/correction.pdf")

    return t_popt, t_pcov


def plot_comp(M):

    filenames = open("/Users/williamquinn/Desktop/set_4/S95_A25/filenames.txt", "r")
    fl = filenames.readlines()
    files = []
    for index, line in enumerate(fl):
        files.append(line.strip())
    path = "/Users/williamquinn/Desktop/set_4/S95_A25/"

    aan_ch0 = []
    aan_ch1 = []
    dates_ch0 = []
    dates_ch1 = []
    for i_file in files:
        date = i_file.split("_")[0]
        try:
            root_file = ROOT.TFile(path + i_file, "READ")
            num_hist = root_file.Get(date + "_GAO607_he_apulse_num_1400V")

            contents = []
            weights = []
            for i_bin in range(num_hist.GetNbinsX() + 1):
                contents.append(i_bin - 1)
                weights.append(num_hist.GetBinContent(i_bin) / num_hist.GetEntries())
            average = np.average(contents, weights=weights)
            aan_ch0.append(average)
            # print(num_hist.GetMean(), average)
            del num_hist
            del root_file
            dates_ch0.append(int(date))
        except:
            print(0, path + i_file, date)
        try:
            root_file_1 = ROOT.TFile(path + i_file, "READ")
            num_hist_1 = root_file_1.Get(date + "_GAO612_he_apulse_num_1400V")

            contents = []
            weights = []
            for i_bin in range(num_hist_1.GetNbinsX() + 1):
                contents.append(i_bin - 1)
                weights.append(num_hist_1.GetBinContent(i_bin) / num_hist_1.GetEntries())
            average = np.average(contents, weights=weights)
            aan_ch1.append(average)
            # print(num_hist.GetMean(), average)
            del num_hist_1
            del root_file_1
            dates_ch1.append(int(date))
        except:
            print(1, path + i_file, date)

    dates_ch0s = process_date(dates_ch0)
    dates_ch1s = process_date(dates_ch1)
    aan_ch0s = np.array(aan_ch0)
    aan_ch1s = np.array(aan_ch1)

    start_ch0 = np.where(dates_ch0s == 98)[0][0]
    start_ch1 = np.where(dates_ch1s == 98)[0][0]
    x = np.array(dates_ch0s[start_ch0:]) - dates_ch0s[start_ch0]
    y = np.array(aan_ch0s[start_ch0:])

    t_popt, t_pcov = adjust_pars(M)

    y_ = t_popt[0] * y * y + t_popt[1] * y + t_popt[2]

    plt.figure(figsize=figsize, facecolor='white')
    plt.plot(x, y, 'C0s', label='Raw Data', markersize=marker_size/2)
    plt.plot(x, y_, 'C2o', label='Corrected Data', markersize=marker_size/2)
    plt.xlim(x[0], x[-1])
    plt.xlabel('Time /days')
    plt.ylabel('Average After-Pulse Number')
    plt.title('Corrected Average After-Pulse Number Comparison')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/PMT_Project/comp.pdf")


def main():
    df = pd.read_csv("/Users/williamquinn/Desktop/PMT_Project/HeMatrix.txt", index_col=0)
    df = df / 500
    df.loc["i0", "0"] = 1
    M = df.values
    plot_comp(M)
    # M = M.T
    df = df.transpose()

    fig = plt.figure(figsize=figsize)
    # fig = plt.subplots()
    p = plt.pcolor(df, cmap='coolwarm')
    cbar = plt.colorbar(p)

    cbar.set_label("Proportion found", rotation=90)
    plt.xlabel("After-Pulses Inserted")
    plt.ylabel("After-Pulses Found")
    plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
    plt.xticks(np.arange(0.5, len(df.index), 1), df.index)
    plt.tight_layout()
    plt.title("Convolution Inefficiency")
    plt.savefig("/Users/williamquinn/Desktop/PMT_Project/HeMatrix.pdf")


if __name__ == "__main__":
    main()
