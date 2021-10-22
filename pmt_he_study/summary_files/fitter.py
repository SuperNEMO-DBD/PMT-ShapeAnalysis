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
    df = pd.read_csv("/Users/williamquinn/Desktop/HeMatrix.txt",
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

    dates_ch0s = process_date(dates_ch0)[:-2]
    dates_ch1s = process_date(dates_ch1)[:-2]
    aan_ch0s = aan_ch0[:-2]
    aan_ch1s = aan_ch1[:-2]

    start_ch0 = np.where(dates_ch0s == 98)[0][0]
    start_ch1 = np.where(dates_ch1s == 98)[0][0]
    x = np.array(dates_ch0s[start_ch0:]) - dates_ch0s[start_ch0]
    y = np.array(aan_ch0s[start_ch0:])

    return x, y


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


def plot_fit(name, model, x, y, yerr, pars, errs, chi):
    fig1 = plt.figure(figsize=(9, 6), facecolor='white')
    frame1 = fig1.add_axes((.1, .3, .8, .6))
    plt.errorbar(x, y, yerr=yerr, fmt='k.', label='Exposed')
    plt.plot(x, model.func(x, pars), 'g-', linewidth=3, zorder=10, label=r'Model $\chi_R=${:.2f}'.format(chi))
    frame1.set_xticklabels([])  # Remove x-tic labels for the first frame
    plt.grid()
    plt.ylabel("Average apulse num")
    plt.xlim(x[0], x[-1])
    plt.legend(loc='upper left')
    plt.title("Average apulse num vs exposure time: " + name)

    frame2 = fig1.add_axes((.1, .1, .8, .2))
    plt.errorbar(x, (model.func(x, pars) - y) / model.func(x, pars), yerr=yerr / model.func(x, pars),
                 fmt='g.')
    plt.axhline(0, ls='--', color='k')
    plt.grid()
    plt.xlim(x[0], x[-1])
    # plt.ylim(-0.25, 0.25)
    plt.ylabel("(model - data)/model")
    plt.xlabel("Days relative to 11/02/2020")
    plt.tight_layout()
    plt.savefig('/Users/williamquinn/Desktop/' + name + '.pdf')
    plt.close()


def main(name: str):
    t_popt, t_pcov = adjust_pars()

    filenames = open("/Users/williamquinn/Desktop/S95_A25/filenames.txt", "r")
    fl = filenames.readlines()
    files = []
    for index, line in enumerate(fl):
        files.append(line.strip())
    path = "/Users/williamquinn/Desktop/S95_A25/"

    x, y = get_data(path, files)
    y_ = t_popt[0] * y * y + t_popt[1] * y + t_popt[0]

    # y = y_

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
    yerr = np.ones_like(y)*0.05

    graph = ROOT.TGraphErrors(len(y), array('d', [i for i in x]), array('d', [i for i in y]),
                              array('d', [1 for i in range(len(x))]), array('d', [i for i in yerr]))

    fit = ROOT.TF1("func", model, float(x[0]), float(x[-1]), len(guess))
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

    plot_fit(name, model, x, y, yerr, pars, errs, chi)

    '''graph.Draw("P")

    graph.SetLineColor(4)
    graph.SetLineWidth(4)
    graph.SetMarkerColor(4)
    graph.SetMarkerSize(0.5)
    graph.SetMarkerStyle(20)
    graph.SetTitle(name)
    graph.GetXaxis().SetTitle("Days since 11/02/2020")
    graph.GetXaxis().SetRangeUser(x[0], x[-1])
    graph.GetYaxis().SetTitle("Average After-pulse number")
    graph.GetYaxis().SetRangeUser(0.0, 2.5)
    graph.Draw("AP")

    c1.Update()
    c1.Modified()
    c1.SaveAs('~/Desktop/' + name + '.pdf')'''

    return


if __name__ == "__main__":
    names = ['Model', 'Model_Eff', 'Model_0', 'Model_Eff_0']
    for name in names:
        main(name)
