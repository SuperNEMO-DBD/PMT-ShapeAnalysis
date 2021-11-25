import sys

sys.path.insert(1, '../..')
import ROOT
import random
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from functions.other_functions import process_date, process_exposure, chi2, linear
import pandas as pd
from array import array
from pmt_he_study.format_plot import *


class Model:
    def __call__(self, x, pars):
        p = 101325 / 10
        p0 = pars[0]
        p1 = pars[1]
        L = pars[2]
        t = x[0] * 3600 * 24

        temp = 0
        for n in range(1, 50):
            temp += ((-1) ** n / n ** 2) * (1 - np.exp(-(n ** 2) * (np.pi ** 2) * t / (L * 6)))
        f2 = (12 / np.pi ** 2) * L * temp
        y = p0 * p * (t + f2) + p1
        return y

    def func(self, x, pars):
        p = 101325 / 10
        p0 = pars[0]
        p1 = pars[1]
        L = pars[2]
        t = x * 3600 * 24

        temp = 0
        for n in range(1, 50):
            temp += ((-1) ** n / n ** 2) * (1 - np.exp(-(n ** 2) * (np.pi ** 2) * t / (L * 6)))
        f2 = (12 / np.pi ** 2) * L * temp
        y = p0 * p * (t + f2) + p1
        return y


class Model_Eff:
    def __call__(self, x, pars):
        p = 101325 / 10
        p0 = pars[0]
        p1 = pars[1]
        L = pars[2]
        A = pars[3]
        B = pars[4]
        t = x[0] * 3600 * 24

        temp = 0
        for n in range(1, 50):
            temp += ((-1) ** n / n ** 2) * (1 - np.exp(-(n ** 2) * (np.pi ** 2) * t / (L * 6)))
        f2 = (12 / np.pi ** 2) * L * temp
        y = (1 - 1/( 1 + np.exp(-1*(t-A)/B)))*(p0 * p * (t + f2) + p1)
        return y

    def func(self, x, pars):
        p = 101325 / 10
        p0 = pars[0]
        p1 = pars[1]
        L = pars[2]
        A = pars[3]
        B = pars[4]
        t = x * 3600 * 24

        temp = 0
        for n in range(1, 50):
            temp += ((-1) ** n / n ** 2) * (1 - np.exp(-(n ** 2) * (np.pi ** 2) * t / (L * 6)))
        f2 = (12 / np.pi ** 2) * L * temp
        y = (1 - 1 / (1 + np.exp(-1 * (t - A) / B))) * (p0 * p * (t + f2) + p1)
        return y


class Model_0:
    def __call__(self, x, pars):
        pe_0 = 101325 / 100
        pe_1 = 101325 / 10
        p0 = pars[0]
        p1 = pars[1]
        L = pars[2]

        t_0 = (x[0] + 98) * 3600 * 24
        t_1 = t_0 - (98 * 3600 * 24)
        temp_0 = 0
        temp_1 = 0
        for n in range(1, 50):
            temp_0 += ((-1) ** n / n ** 2) * (1 - np.exp(-(n ** 2) * (np.pi ** 2) * t_0 / (L * 6)))
            temp_1 += ((-1) ** n / n ** 2) * (1 - np.exp(-(n ** 2) * (np.pi ** 2) * t_1 / (L * 6)))
        f_0 = (12 / np.pi ** 2) * L * temp_0
        f_1 = (12 / np.pi ** 2) * L * temp_1

        y = p0 * (pe_0 * (t_0 + f_0) + pe_1 * (t_1 + f_1)) + p1
        return y

    def func(self, x, pars):
        pe_0 = 101325 / 100
        pe_1 = 101325 / 10
        p0 = pars[0]
        p1 = pars[1]
        L = pars[2]

        t_0 = (x + 98) * 3600 * 24
        t_1 = t_0 - (98 * 3600 * 24)
        temp_0 = 0
        temp_1 = 0
        for n in range(1, 50):
            temp_0 += ((-1) ** n / n ** 2) * (1 - np.exp(-(n ** 2) * (np.pi ** 2) * t_0 / (L * 6)))
            temp_1 += ((-1) ** n / n ** 2) * (1 - np.exp(-(n ** 2) * (np.pi ** 2) * t_1 / (L * 6)))
        f_0 = (12 / np.pi ** 2) * L * temp_0
        f_1 = (12 / np.pi ** 2) * L * temp_1

        y = p0 * (pe_0 * (t_0 + f_0) + pe_1 * (t_1 + f_1)) + p1
        return y


class Model_Eff_0:
    def __call__(self, x, pars):
        pe_0 = 101325 / 100
        pe_1 = 101325 / 10
        p0 = pars[0]
        p1 = pars[1]
        L = pars[2]
        A = pars[3]
        B = pars[4]

        t_0 = (x[0] + 98) * 3600 * 24
        t_1 = t_0 - (98 * 3600 * 24)
        temp_0 = 0
        temp_1 = 0
        for n in range(1, 50):
            temp_0 += ((-1) ** n / n ** 2) * (1 - np.exp(-(n ** 2) * (np.pi ** 2) * t_0 / (L * 6)))
            temp_1 += ((-1) ** n / n ** 2) * (1 - np.exp(-(n ** 2) * (np.pi ** 2) * t_1 / (L * 6)))
        f_0 = (12 / np.pi ** 2) * L * temp_0
        f_1 = (12 / np.pi ** 2) * L * temp_1

        y = (1 - 1/( 1 + np.exp(-1*(t_1-A)/B)))*(p0 * (pe_0 * (t_0 + f_0) + pe_1 * (t_1 + f_1)) + p1)
        return y

    def func(self, x, pars):
        pe_0 = 101325 / 100
        pe_1 = 101325 / 10
        p0 = pars[0]
        p1 = pars[1]
        L = pars[2]
        A = pars[3]
        B = pars[4]

        t_0 = (x + 98) * 3600 * 24
        t_1 = t_0 - (98 * 3600 * 24)
        temp_0 = 0
        temp_1 = 0
        for n in range(1, 50):
            temp_0 += ((-1) ** n / n ** 2) * (1 - np.exp(-(n ** 2) * (np.pi ** 2) * t_0 / (L * 6)))
            temp_1 += ((-1) ** n / n ** 2) * (1 - np.exp(-(n ** 2) * (np.pi ** 2) * t_1 / (L * 6)))
        f_0 = (12 / np.pi ** 2) * L * temp_0
        f_1 = (12 / np.pi ** 2) * L * temp_1

        y = (1 - 1 / (1 + np.exp(-1 * (t_1 - A) / B))) * (p0 * (pe_0 * (t_0 + f_0) + pe_1 * (t_1 + f_1)) + p1)
        return y


def my_matrix():
    df = pd.read_csv("/Users/williamquinn/Desktop/PMT_Project/HeMatrix.txt",
                     index_col=0)
    df = df / 500
    df.loc["i0", "0"] = 1
    M = df.values

    return M


def get_data(path, files):
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
    aan_ch0s = aan_ch0
    aan_ch1s = aan_ch1

    start_ch0 = np.where(dates_ch0s == 98)[0][0]
    start_ch1 = np.where(dates_ch1s == 98)[0][0]

    x0 = np.array(dates_ch0s[start_ch0:]) - dates_ch0s[start_ch0]
    y0 = np.array(aan_ch0s[start_ch0:])
    x1 = np.array(dates_ch1s[start_ch1:]) - dates_ch1s[start_ch1]
    y1 = np.array(aan_ch1s[start_ch1:])

    return x0, y0, x1, y1


def test_model(x, q, m, c):
    return q * x * x + m * x + c


def adjust_pars():
    averages = []
    stds = []
    M = my_matrix()
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

    means = np.linspace(0, 18, 19 * 5)
    f = []
    for i_mean in means:
        probs = []
        for n in nums:
            prob = np.exp(-1 * i_mean) * (i_mean ** n) / np.math.factorial(n) * averages[n]
            probs.append(prob)
        f.append(np.sum(probs))

    f = np.array(f)
    my_list = np.where(means < 5)[0]

    t_popt, t_pcov = curve_fit(f=test_model, xdata=f[my_list], ydata=means[my_list])
    return t_popt, t_pcov


def plot_fit(name, model, x0, y0, yerr0, x1, y1, yerr1, pars, errs, chi):
    fig1 = plt.figure(figsize=figsize, facecolor='white')
    frame1 = fig1.add_axes((.13, .35, .8, .6))
    frame1.set_xticklabels([])

    plt.errorbar(x0, y0, yerr=yerr0, fmt='C0.', label='Exposed', zorder=0, markersize=1, linewidth=line_width)
    plt.errorbar(x1, y1, yerr=yerr1, fmt='C1.', label='Control', zorder=0, markersize=1, linewidth=line_width)
    plt.plot(x0, model.func(x0, pars), 'C4-', linewidth=line_width*2, zorder=10, label=r'Model')
    plt.ylabel("Average After-Pulse Number")
    plt.xlim(x0[0], x0[-1])

    handles, labels = plt.gca().get_legend_handles_labels()

    patch = patches.Patch(color='white', label=r'$P_0 =$ {:.2e} ± {:.0e}'.format(pars[0], errs[0]))
    patch_1 = patches.Patch(color='white', label=r'$P_1 =$ {:.2f} ± {:.2f}'.format(pars[1], errs[1]))
    patch_2 = patches.Patch(color='white', label=r'$L =$ {:.0f} ± {:.0f}'.format(pars[2]/(3600*24), errs[2]/(3600*24)))
    patch_3 = patches.Patch(color='white', label=r'Model $\chi_R=${:.2f}'.format(chi))

    handles.extend([patch, patch_1, patch_2, patch_3])
    plt.legend(handles=handles, loc='best')
    plt.title("Average Apulse Number vs Exposure Time")

    frame2 = fig1.add_axes((.13, .15, .8, .2))
    plt.errorbar(x0, (model.func(x0, pars) - y0) / model.func(x0, pars), yerr=yerr0 / model.func(x0, pars),
                 fmt='C2.', markersize=1, zorder=0, linewidth=line_width)
    plt.axhline(0, ls='--', color='k', zorder=10, linewidth=line_width)
    plt.xlim(x0[0], x0[-1])
    # plt.ylim(-0.25, 0.25)
    plt.ylabel("(model - data)/model", fontsize=4)
    plt.xlabel("Days relative to 11/02/2020")
    plt.tight_layout()
    plt.savefig('/Users/williamquinn/Desktop/' + name + '.pdf')
    plt.close()


def main(name: str):
    t_popt, t_pcov = adjust_pars()

    filenames = open("/Users/williamquinn/Desktop/set_4/S95_A25/filenames.txt", "r")
    fl = filenames.readlines()
    files = []
    for index, line in enumerate(fl):
        files.append(line.strip())
    path = "/Users/williamquinn/Desktop/set_4/S95_A25/"

    x0, y0, x1, y1 = get_data(path, files)
    y_0 = t_popt[0] * y0 * y0 + t_popt[1] * y0 + t_popt[0]
    y_1 = t_popt[0] * y1 * y1 + t_popt[1] * y1 + t_popt[0]

    '''c1 = ROOT.TCanvas("c1", "c1", 200, 10, 600, 400)
    c1.SetGrid()'''

    if name == 'Model':
        guess = array('d', [1.52157116e-11, 6.84547311e-02, 1.13069872e+07])
        names = ["p0", "p1", "L"]
        model = Model()
    elif name == 'Model_Eff':
        guess = array('d', [2.73836056e-11, 1.06234990e-01, 1.68440129e+07, 3.18837265e+07, 5.73493526e+06])
        names = ["p0", "p1", "L", "A", "B"]
        model = Model_Eff()
    elif name == 'Model_0':
        guess = array('d', [1.52157116e-11, 6.84547311e-02, 1.13069872e+07])
        names = ["p0", "p1", "L"]
        model = Model_0()
    elif name == 'Model_Eff_0':
        guess = array('d', [2.73836056e-11, 1.06234990e-01, 1.68440129e+07, 3.18837265e+07, 5.73493526e+06])
        names = ["p0", "p1", "L", "A", "B"]
        model = Model_Eff_0()
    else:
        return
    yerr0 = np.ones_like(y0) * 0.05
    yerr1 = np.ones_like(y1) * 0.05
    yerr_0 = np.ones_like(y_0) * 0.05
    yerr_1 = np.ones_like(y_1) * 0.05

    graph = ROOT.TGraphErrors(len(y0), array('d', [i for i in x0]), array('d', [i for i in y0]),
                              array('d', [1 for i in range(len(x0))]), array('d', [i for i in yerr0]))

    fit = ROOT.TF1("func", model, float(x0[0]), float(x0[-1]), len(guess))
    for i in range(len(guess)):
        fit.SetParameter(i, guess[i])
        fit.SetParName(i, names[i])
        # fit.SetParLimits(i, guess[i] - abs(guess[i]), guess[i] + abs(guess[i]))

    graph.Fit("func", "R")
    pars = []
    errs = []
    chi = fit.GetChisquare()/fit.GetNDF()

    for i in range(len(guess)):
        pars.append(fit.GetParameter(i))
        errs.append(fit.GetParError(i))
    del graph
    del fit

    plot_fit('raw_aan_data', model, x0, y0, yerr0, x1, y1, yerr1, pars, errs, chi)

    graph = ROOT.TGraphErrors(len(y_0), array('d', [i for i in x0]), array('d', [i for i in y_0]),
                              array('d', [1 for i in range(len(x0))]), array('d', [i for i in yerr_0]))

    fit = ROOT.TF1("func", model, float(x0[0]), float(x0[-1]), len(guess))
    for i in range(len(guess)):
        fit.SetParameter(i, guess[i])
        fit.SetParName(i, names[i])
        # fit.SetParLimits(i, guess[i] - abs(guess[i]), guess[i] + abs(guess[i]))

    graph.Fit("func", "R")
    pars = []
    errs = []
    chi = fit.GetChisquare() / fit.GetNDF()

    for i in range(len(guess)):
        pars.append(fit.GetParameter(i))
        errs.append(fit.GetParError(i))
    del graph
    del fit

    plot_fit('adjusted_aan_data', model, x0, y_0, yerr_0, x1, y_1, yerr_1, pars, errs, chi)

    return


if __name__ == "__main__":
    main('Model')
