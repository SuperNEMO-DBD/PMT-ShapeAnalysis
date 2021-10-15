import ROOT
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import fnmatch
from datetime import date
import pandas as pd
from scipy.special import factorial
from scipy.optimize import curve_fit
from datetime import date
from scipy.stats import chisquare
import tqdm


def directory_list(folder_path):
    files = []
    filenames = open(folder_path + "/filenames.txt", "r")
    fl = filenames.readlines()
    for index, line in enumerate(fl):
        if len(line.split(".")) > 0:
            filename = line.strip()
            files.append(filename)
        else:
            continue
    return files
    '''if os.path.exists(folder_path):
        for filename in os.listdir(folder_path):
            if fnmatch.fnmatch(filename, '*1400*'):
                files.append(filename)
        return files
    else:
        print("Folder not found")
        sys.exit()'''


def gaussian(x, mu, sd, A):
    exponent = -0.5 * np.square((x - mu) / sd)
    return A * np.exp(exponent)


def gaus_mod(x, mu, sd, A):
    diff = x - mu
    exponent = np.exp(-0.5 * (diff / sd) ** 2)
    dist = A * exponent
    return dist_smear(dist)


def day_extractor(dates):
    start = str(np.min(dates))
    start_date = date(int("20" + start[0:2]), int(start[2:4]), int(start[4:6]))
    days = []

    for i_date in dates:
        i_date = str(i_date)
        current_date = date(int("20" + i_date[0:2]), int(i_date[2:4]), int(i_date[4:6]))
        difference = current_date - start_date
        days.append(difference.days)

    return days


def get_pulses(path, filename):
    # creating the file path
    full_path = path + filename

    # Opening the root file
    root_file = ROOT.TFile(full_path, "READ")
    root_file.cd()

    hist = root_file.Get("110011_GAO607_he_apulse_times_1400V")
    total = hist.GetEntries()

    return total


def get_bins_values(path):
    filename = path.split("/")[-1]
    date = filename.split("_")[0]

    root_file = ROOT.TFile(path, "READ")
    root_file.cd()

    hist = root_file.Get(date + "_GAO607_apulse_num_1400V")
    mean = hist.GetMean()

    bins = []
    freq = []
    values = []

    for i in range(1, hist.GetNbinsX()):
        bins.append(int(i - 1))
        freq.append(int(hist.GetBinContent(i)))

    bins = np.array(bins)
    freq = np.array(freq)
    normed = freq / sum(freq)

    return bins, normed, mean, freq


def dist_smear(start_values):
    '''df = pd.read_csv("/unix/nemo4/PMT_He_Study_nemo4/PMT-ShapeAnalysis/pmt_he_study/waveform_sims/Matrix.txt",
                     index_col=0)'''
    df = pd.read_csv("/unix/nemo4/PMT_He_Study_nemo4/PMT-ShapeAnalysis/pmt_he_study/waveform_sims/HeMatrix.txt",
                     index_col=0)
    df = df / 500
    df.loc["i0", "0"] = 1
    M = df.values
    transformed_values = M.T @ start_values
    return transformed_values


def mean_calc(bins, dist):
    dot_prod = np.dot(bins, dist)
    return dot_prod / sum(dist)


def mychi(O, Exp, Err):
    chi = 0
    parameters = 1
    bins = 0

    for i in range(len(O)):
        if O[i] == 0:
            None
        else:
            bins = bins + 1
            chi = chi + (O[i] - Exp[i]) ** 2 / Err[i] ** 2

    return chi / (bins - parameters)


'''
def func_form(t,A,B,t_onset):
  exp = np.exp(t-t_onset)
  denom = 1 + B*exp
  return A/denom
'''


def func_form(t, A, B, C):
    exp = np.exp(-C * t)
    return A * (1 - 1 / (1 + B * exp))


def new_func_form(t, A, B, C):
    exp = np.exp(-C * (t - B))
    return A * (1 - 1 / (1 + exp))


def main(path, name):
    files = directory_list(path)
    # del files[0]
    dates = []
    ratio = []

    for i_file in tqdm.tqdm(files):
        split = i_file.split("_")
        i_date = split[0]
        if "." in i_date:
            None
        else:
            if int(i_date) < 200211 or i_date == "200317":
                None
            else:

                bins, real_values, real_mean, freq = get_bins_values(path + i_file)
                bins = bins + 0.5  # the bin centre is the one that needs to be used

                freq = [0.000001 if x == 0 else x for x in freq]
                err = np.sqrt(freq) / sum(freq)

                try:
                    p_guess = [0, 1, 1]
                    p_bounds = [[0, 0, 0], [20, 20, 10]]
                    popt, pcov = curve_fit(gaus_mod, bins, real_values, sigma=err, maxfev=2000,
                                           p0=p_guess, bounds=p_bounds)
                except:
                    print("{} Not fitted".format(i_date))
                    continue

                new_real_mean = popt[0]
                sd = popt[1]
                A = popt[2]

                start_values = gaussian(bins, new_real_mean, sd, A)
                transformed_values = dist_smear(start_values)

                dist_mean = mean_calc(bins, start_values)
                smear_mean = mean_calc(bins, transformed_values)
                # print(real_mean, new_real_mean, dist_mean, smear_mean)

                ratio.append(smear_mean / dist_mean)
                dates.append(int(i_date))

                # print("{} Done".format(i_date), "Ratio", smear_mean / dist_mean)

                x = np.linspace(0, 19, 100)

                '''fig = plt.figure(figsize=(9, 6), facecolor='white')
                plt.bar(bins, real_values, width=1, color='blue', label='data')
                plt.plot(bins, start_values, "r.", label='start values')
                plt.plot(x, gaussian(x, new_real_mean, sd, A), "r--")
                plt.plot(bins, transformed_values, "g.", label='smeared values')
                plt.grid()
                plt.xlabel("apulse number")
                plt.ylabel("normalised counts")
                plt.legend(loc='upper right')
                plt.savefig("plots/{}_apnum_fit.pdf".format(i_date))
                plt.close()'''

    st_date = str(np.min(dates))
    days = day_extractor(dates)

    days = np.array(days)
    ratio = np.array(ratio)

    start_date = date(int("20" + st_date[0:2]), int(st_date[2:4]), int(st_date[4:6]))
    start_date = date(int("20" + st_date[0:2]), int(st_date[2:4]), int(st_date[4:6]))

    days_of_ten = date(2020, 2, 11) - start_date

    error = ratio * 0.01

    '''popt, pcov = curve_fit(func_form, days, ratio, maxfev=2000, sigma=error)
    A = popt[0]
    B = popt[1]
    C = popt[2]'''

    popt, pcov = curve_fit(new_func_form, days, ratio, maxfev=2000, sigma=error,
                           p0=[0.9, ],
                           bounds=[[], []])
    A = popt[0]
    B = popt[1]
    C = popt[2]
    D = popt[3]

    # A = 0.97
    # B = 11645
    # C = 0.01

    for i_popt in popt:
        print(i_popt)

    print(pcov)

    '''for g in pcov:
        print(g)'''

    days = np.sort(days)

    # fitting = func_form(days, A, B, C)
    fitting = new_func_form(days, *popt)
    chi_value = mychi(fitting, ratio, error) / 3
    residuals = ratio - fitting

    fig, (ax0, ax1) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]})

    ax0.set_title("Efficiency factor from He simulation from 10% flow")
    ax0.scatter(days, ratio, s=1, color='b', label="Calculated ratios")
    ax0.errorbar(days, ratio, yerr=error, color='b', fmt='.', elinewidth=1)
    ax0.plot(days, fitting, alpha=1, color='r', label="Fitting")
    # ax0.axvline(days_of_ten.days, label="10% Helium", color = "orange")

    # ax0.annotate(r"$\chi^2_{R}$ :"  + str(round(chi_value,2)), xy=(0.01,0.6), xycoords="axes fraction", horizontalalignment ="left")
    # ax0.annotate("A: "  + str(round(A,2)), xy=(0.01,0.5), xycoords="axes fraction",horizontalalignment ="left")
    # ax0.annotate("B: "  + str(round(B,2)), xy=(0.01,0.4), xycoords="axes fraction",horizontalalignment ="left")
    # ax0.annotate("C: "  + str(round(C,2)), xy=(0.01,0.3), xycoords="axes fraction",horizontalalignment ="left")

    ax0.set_ylabel("Efficiency factor")
    ax0.legend(loc="lower left")
    ax0.grid()

    ax1.scatter(days, residuals, color="k", s=1)
    ax1.set_xlabel("Days since " + st_date[4:6] + "/" + st_date[2:4] + "/" + st_date[0:2])
    ax1.set_ylabel("Residuals")
    # ax1.hlines(0,xmin=0,xmax=, color="y", linestyles="dashed")
    ax1.axhline(0, color="y", ls="--")
    ax1.grid()

    fig.tight_layout()
    fig.savefig(name + ".pdf")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Must enter folder path and output image name")
        print("Example 'python SmearRatioWithExposure.py nemo/rootfiles/ final_image'")
    else:
        folder_path = sys.argv[1]
        image_name = sys.argv[2]
        main(folder_path, image_name)
