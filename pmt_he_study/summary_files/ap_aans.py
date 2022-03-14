import sys

sys.path.insert(1, '../..')

import ROOT
import tqdm
import numpy as np
import pandas as pd
from array import array
from scipy.optimize import curve_fit
from pmt_he_study.format_plot import *
from functions.other_functions import *
from src.PMT_Classes import *


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


def my_matrix():
    df = pd.read_csv("/Users/williamquinn/Desktop/PMT_Project/HeMatrix.txt",
                     index_col=0)
    df = df / 500
    df.loc["i0", "0"] = 1
    M = df.values

    return M


def test_model(x, q, m, c):
    return q * x * x + m * x + c


def root_fit(x, y, yerr, model):
    # guess = array('d', [1.52157116e-11, 6.84547311e-02, 1.13069872e+07])
    guess = array('d', [(y[-1] - y[-50])/((x[-1] - x[-50])*3600*24) /10000, y[0], 300*3600*24])
    names = ["p0", "p1", "L"]
    graph = ROOT.TGraphErrors(len(y), array('d', [i for i in x]), array('d', [i for i in y]),
                              array('d', [1 for i in range(len(x))]), array('d', [i for i in yerr]))

    fit = ROOT.TF1("func", model, float(x[0]), float(x[-1]), len(guess))
    for i in range(len(guess)):
        fit.SetParameter(i, guess[i])
        fit.SetParName(i, names[i])
        fit.SetParLimits(i, guess[i] - abs(guess[i])/2, guess[i] + abs(guess[i])/2)

    graph.Fit("func", "0Q")
    pars = []
    errs = []
    chi = fit.GetChisquare() / fit.GetNDF()

    for i in range(len(guess)):
        pars.append(fit.GetParameter(i))
        errs.append(fit.GetParError(i))
    return pars, errs, chi


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


def plot_aan(dates, aans, aans_err, name, model):
    # model = Model()
    t_popt, t_pcov = adjust_pars()

    date_0 = process_date(dates[0])
    date_1 = process_date(dates[1])
    try:
        start = np.where(date_0 == 0)[0][0]
    except:
        start = np.where(date_0 == 1)[0][0]
    mid = np.where(date_0 == 98)[0][0]

    aan_0 = np.array(aans[0])
    aan_1 = np.array(aans[1])

    aan_err_0 = np.array(aans_err[0])
    aan_err_1 = np.array(aans_err[1])

    if name == "he":
        aan_0 = t_popt[0] * aan_0 * aan_0 + t_popt[1] * aan_0 + t_popt[2]
        aan_err_0 = aan_err_0*2
        aan_1 = t_popt[0] * aan_1 * aan_1 + t_popt[1] * aan_1 + t_popt[2]
        aan_err_1 = aan_err_1*2

    x = np.array(date_0[mid + 1:]) - date_0[mid + 1:][0]
    y = np.array(aan_0[mid + 1:])
    yerr = np.array(aan_err_0[mid + 1:])
    pars, errs, chi = root_fit(x, y, yerr, model)

    fig1 = plt.figure(num=None, figsize=(5, 5), dpi=80, facecolor='w', edgecolor='k')
    frame1 = fig1.add_axes((.12, .32, .8, .6))
    frame1.set_xticklabels([])
    plt.errorbar(date_0[:start + 1], aan_0[:start + 1], zorder=0,
                 yerr=aan_err_0[:start + 1], fmt="C0s", label="Atmospheric He", markersize=1)
    plt.errorbar(date_0[start + 1:mid + 1], aan_0[start + 1:mid + 1], zorder=0,
                 yerr=aan_err_0[start + 1:mid + 1], fmt="C1s", label="1% He", markersize=1)
    plt.errorbar(date_0[mid + 1:], aan_0[mid + 1:], zorder=0,
                 yerr=aan_err_0[mid + 1:], fmt="C2s", label="10% He", markersize=1)
    plt.errorbar(date_1, aan_1, zorder=0,
                 yerr=aan_err_1, fmt="C3o", label="Control", markersize=1)
    plt.plot(date_0[mid + 1:], model.func(x, pars), 'k-', label='Model', zorder=10)

    handles, labels = plt.gca().get_legend_handles_labels()
    patch = patches.Patch(color='white', label=r'$P_0 =$ {:.3e} ± {:.0e}'.format(pars[0], errs[0]))
    patch_1 = patches.Patch(color='white', label=r'$P_1 =$ {:.3e} ± {:.0e}'.format(pars[1], errs[1]))
    patch_2 = patches.Patch(color='white', label=r'$L =$ {:.0f} ± {:.0f}'.format(pars[2] / (3600 * 24), errs[2] / (3600 * 24)))
    patch_3 = patches.Patch(color='white', label=r'$\chi^2_R =$ {:.2f}'.format(chi))
    handles.extend([patch, patch_1, patch_2, patch_3])

    plt.ylabel("Average After-Pulse Number")
    plt.title("Average After-Pulse Number vs Exposure")
    plt.xlim(-30, 420)
    plt.legend(handles=handles, loc='upper left')

    frame2 = fig1.add_axes((.12, .1, .8, .2))
    plt.xlabel("Days from 1% He Onset")
    plt.axhline(0, ls='--', color='black')
    plt.ylabel("(model-data)/model")
    plt.xlim(-30, 420)
    plt.errorbar(date_0[mid + 1:],
                 (model.func(x, pars) - aan_0[mid + 1:]) / model.func(x, pars),
                 yerr=aan_err_0[mid + 1:]/model.func(date_0[mid + 1:], pars), fmt="k.")
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/aan_vs_time_" + name + ".pdf")
    plt.close()


def main():
    topology = [2, 1]
    pmt_array = PMT_Array(topology, "summary")
    pmt_array.set_pmt_id("GAO607", 0)
    pmt_array.set_pmt_id("GAO612", 1)

    filenames_txt = "/Users/williamquinn/Desktop/data/1400V/filenames.txt"
    try:
        print(">>> Reading data from file: {}".format(filenames_txt))
        date_file = open(filenames_txt, 'r')
    except FileNotFoundError as fnf_error:
        print(fnf_error)
        raise Exception("Error opening data file {}".format(filenames_txt))

    filenames = np.loadtxt(filenames_txt, delimiter=',', dtype={
        'names': ['filename'],
        'formats': ['S100']}, unpack=True)

    dates = [[], []]
    nums = [[], []]
    he_nums = [[], []]
    nums_err = [[], []]
    he_nums_err = [[], []]

    for i_file in tqdm.tqdm(range(filenames.size)):
        filename = filenames[i_file][0].decode("utf-8")
        date = filename.split("_")[0]
        voltage = int(filename.split("_")[1].split("A")[1])

        file = ROOT.TFile("/Users/williamquinn/Desktop/data/1400V/" + filename, "READ")
        file.cd()

        for i_om in range(2):
            num_hist = file.Get(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() +
                                      "_apulse_num_" + str(voltage) + "V")
            he_num_hist = file.Get(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() +
                                      "_he_apulse_num_" + str(voltage) + "V")

            try:
                num_hist.GetEntries()
                he_num_hist.GetEntries()
            except:
                continue

            dates[i_om].append(int(date))

            nums[i_om].append(num_hist.GetMean())
            he_nums[i_om].append(he_num_hist.GetMean())

            nums_err[i_om].append(num_hist.GetMeanError())
            he_nums_err[i_om].append(he_num_hist.GetMeanError())

            del num_hist
            del he_num_hist

    model = Model_0()
    plot_aan(dates, nums, nums_err, "", model)
    plot_aan(dates, he_nums, he_nums_err, "he", model)


if __name__ == "__main__":
    main()
