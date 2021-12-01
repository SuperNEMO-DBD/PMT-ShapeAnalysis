import sys

sys.path.insert(1, '../..')

import ROOT
import tqdm
import numpy as np
from array import array
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


def plot_ap_charge(dates, ap_charges, ap_charges_err, name):
    date_0 = process_date(dates[0])
    date_1 = process_date(dates[1])
    try:
        start = np.where(date_0 == 0)[0][0]
    except:
        start = np.where(date_0 == 1)[0][0]
    mid = np.where(date_0 == 98)[0][0]

    ap_charge_0 = np.array(ap_charges[0])
    ap_charge_1 = np.array(ap_charges[1])

    ap_charge_err_0 = np.array(ap_charges_err[0])
    ap_charge_err_1 = np.array(ap_charges_err[1])

    x = np.array(date_0[mid + 1:]) - date_0[mid + 1:][0]
    y = np.array(ap_charge_0[mid + 1:])
    yerr = np.array(np.array(ap_charge_err_0[mid + 1:]))
    pars, errs, chi = root_fit(x, y, yerr)

    fig1 = plt.figure(num=None, figsize=(5,5), dpi=80, facecolor='w', edgecolor='k')
    frame1 = fig1.add_axes((.1, .3, .8, .6))
    frame1.set_xticklabels([])
    plt.errorbar(date_0[:start + 1], ap_charge_0[:start + 1], zorder=0,
                 yerr=ap_charge_err_0[:start + 1], fmt="C0s", label="Atmospheric He", markersize=1)
    plt.errorbar(date_0[start + 1:mid + 1], ap_charge_0[start + 1:mid + 1], zorder=0,
                 yerr=ap_charge_err_0[start + 1:mid + 1], fmt="C1s", label="1% He", markersize=1)
    plt.errorbar(date_0[mid + 1:], ap_charge_0[mid + 1:], zorder=0,
                 yerr=ap_charge_err_0[mid + 1:], fmt="C2s", label="10% He", markersize=1)
    plt.errorbar(date_1, ap_charge_1, zorder=0,
                 yerr=ap_charge_err_1, fmt="C3o", label="Control", markersize=1)
    plt.plot(date_0[mid + 1:], Model().func(x, pars), 'k-', label='Model', zorder=10)
    # plt.axvline(date_0[start], 0, 100, ls='--', color='k')
    # plt.axvline(date_0[mid], 0, 100, ls='--', color='k')

    handles, labels = plt.gca().get_legend_handles_labels()
    patch = patches.Patch(color='white', label=r'$P_0 =$ {:.4e} ± {:.0e}'.format(pars[0], errs[0]))
    patch_1 = patches.Patch(color='white', label=r'$P_1 =$ {:.3f} ± {:.3f}'.format(pars[1], errs[1]))
    patch_2 = patches.Patch(color='white', label=r'$L =$ {:.0f} ± {:.0f}'.format(pars[2] / (3600 * 24), errs[2] / (3600 * 24)))
    patch_3 = patches.Patch(color='white', label=r'$\chi^2_R =$ {:.2f}'.format(chi))
    handles.extend([patch, patch_1, patch_2, patch_3])

    plt.ylabel("Average Charge /pC")
    plt.title("Average After-pulse Region Charge")
    plt.xlim(-30, 420)
    plt.legend(handles=handles, loc='upper left')

    frame2 = fig1.add_axes((.1, .1, .8, .2))
    plt.xlabel("Exposure Days Relative to 06/11/2019")
    plt.axhline(0, ls='--', color='black')
    plt.ylabel("(model-data)/model")
    plt.xlim(-30, 420)
    plt.errorbar(date_0[mid + 1:],
                 (Model().func(x, pars) - ap_charge_0[mid + 1:]) / Model().func(x, pars),
                 yerr=ap_charge_err_0[mid + 1:]/Model().func(date_0[mid + 1:], pars), fmt="k.")
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/PMT_Project/ap_charge_vs_time_" + name + ".pdf")
    plt.close()


def plot_aapc_vs_charge(dates, ratios, ratios_err, name: str):
    date_0 = process_date(dates[0])
    date_1 = process_date(dates[1])
    try:
        start = np.where(date_0 == 0)[0][0]
    except:
        start = np.where(date_0 == 1)[0][0]
    mid = np.where(date_0 == 98)[0][0]

    ratio_0 = np.array(ratios[0])
    ratio_1 = np.array(ratios[1])

    ratio_err_0 = np.array(ratios_err[0])
    ratio_err_1 = np.array(ratios_err[1])

    x = np.array(date_0[mid + 1:]) - date_0[mid + 1:][0]
    y = np.array(ratio_0[mid + 1:])
    yerr = np.array(ratio_err_0[mid + 1:])
    pars, errs, chi = root_fit(x, y, yerr)

    fig1 = plt.figure(num=None, figsize=(5,5), dpi=80, facecolor='w', edgecolor='k')
    frame1 = fig1.add_axes((.1, .3, .8, .6))
    frame1.set_xticklabels([])
    plt.errorbar(date_0[:start + 1], ratio_0[:start + 1], zorder=0,
                 yerr=ratio_err_0[:start + 1], fmt="C0s", label="Atmospheric He", markersize=1)
    plt.errorbar(date_0[start + 1:mid + 1], ratio_0[start + 1:mid + 1], zorder=0,
                 yerr=ratio_err_0[start + 1:mid + 1], fmt="C1s", label="1% He", markersize=1)
    plt.errorbar(date_0[mid + 1:], ratio_0[mid + 1:], zorder=0,
                 yerr=ratio_err_0[mid + 1:], fmt="C2s", label="10% He", markersize=1)
    plt.errorbar(date_1, ratio_1, zorder=0,
                 yerr=ratio_err_1, fmt="C3o", label="Control", markersize=1)
    plt.plot(date_0[mid + 1:], Model().func(x, pars), 'k-', label='Model', zorder=10)
    # plt.axvline(date_0[start], 0, 100, ls='--', color='k')
    # plt.axvline(date_0[mid], 0, 100, ls='--', color='k')

    handles, labels = plt.gca().get_legend_handles_labels()
    patch = patches.Patch(color='white', label=r'$P_0 =$ {:.4e} ± {:.0e}'.format(pars[0], errs[0]))
    patch_1 = patches.Patch(color='white', label=r'$P_1 =$ {:.4f} ± {:.4f}'.format(pars[1], errs[1]))
    patch_2 = patches.Patch(color='white',
                            label=r'$L =$ {:.0f} ± {:.0f}'.format(pars[2] / (3600 * 24), errs[2] / (3600 * 24)))
    patch_3 = patches.Patch(color='white', label=r'$\chi^2_R =$ {:.2f}'.format(chi))
    handles.extend([patch, patch_1, patch_2, patch_3])

    plt.ylabel("Charge Ratio Average")
    plt.title("After-pulse Region Charge vs Pulse Charge Ratio")
    plt.xlim(-30, 420)
    plt.legend(handles=handles, loc='upper left')

    frame2 = fig1.add_axes((.1, .1, .8, .2))
    plt.xlabel("Exposure Days Relative to 06/11/2019")
    plt.axhline(0, ls='--', color='black')
    plt.ylabel("(model-data)/model")
    plt.xlim(-30, 420)
    plt.errorbar(date_0[mid + 1:],
                 (Model().func(x, pars) - ratio_0[mid + 1:]) / Model().func(x, pars),
                 yerr=ratio_err_0[mid + 1:] / Model().func(date_0[mid + 1:], pars), fmt="k.")

    plt.savefig("/Users/williamquinn/Desktop/PMT_Project/ratio_vs_time_" + name + ".pdf")
    plt.close()


def root_fit(x, y, yerr):
    # guess = array('d', [1.52157116e-11, 6.84547311e-02, 1.13069872e+07])
    guess = array('d', [(y[-1] - y[-50])/((x[-1] - x[-50])*3600*24) /10000, y[0], 300*3600*24])
    names = ["p0", "p1", "L"]
    model = Model()
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


def main():
    topology = [2, 1]
    pmt_array = PMT_Array(topology, "summary")
    pmt_array.set_pmt_id("GAO607", 0)
    pmt_array.set_pmt_id("GAO612", 1)

    filenames_txt = "/Users/williamquinn/Desktop/set_4/S95_A25/filenames.txt"
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

    ap_charge = [[], []]
    ap_ratio = [[], []]
    he_ap_charge = [[], []]
    he_ap_ratio = [[], []]

    ap_charge_err = [[], []]
    ap_ratio_err = [[], []]
    he_ap_charge_err = [[], []]
    he_ap_ratio_err = [[], []]


    for i_file in tqdm.tqdm(range(filenames.size)):
        filename = filenames[i_file][0].decode("utf-8")
        date = filename.split("_")[0]
        voltage = int(filename.split("_")[1].split("A")[1])

        file = ROOT.TFile("/Users/williamquinn/Desktop/set_4/S95_A25/" + filename, "READ")
        file.cd()

        for i_om in range(2):
            ap_charge_hist = file.Get(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() +
                                      "_ap_charge_spectrum_" + str(voltage) + "V")
            he_ap_charge_hist = file.Get(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() +
                                         "_he_ap_charge_spectrum_" + str(voltage) + "V")
            ap_charge_charge_hist = file.Get(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() +
                                             "_ap_charge_charge_spectrum_" + str(voltage) + "V")
            he_ap_charge_charge_hist = file.Get(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() +
                                                "_he_ap_charge_charge_spectrum_" + str(voltage) + "V")

            try:
                ap_charge_hist.GetEntries()
                he_ap_charge_hist.GetEntries()
                ap_charge_charge_hist.GetEntries()
                he_ap_charge_charge_hist.GetEntries()
            except:
                continue

            dates[i_om].append(int(date))

            ap_charge[i_om].append(ap_charge_hist.GetMean())
            he_ap_charge[i_om].append(he_ap_charge_hist.GetMean())
            ap_ratio[i_om].append(ap_charge_charge_hist.GetMean())
            he_ap_ratio[i_om].append(he_ap_charge_charge_hist.GetMean())

            ap_charge_err[i_om].append(ap_charge_hist.GetMeanError())
            he_ap_charge_err[i_om].append(he_ap_charge_hist.GetMeanError())
            ap_ratio_err[i_om].append(ap_charge_charge_hist.GetMeanError())
            he_ap_ratio_err[i_om].append(he_ap_charge_charge_hist.GetMeanError())

            del ap_charge_hist
            del he_ap_charge_hist
            del ap_charge_charge_hist
            del he_ap_charge_charge_hist

    plot_ap_charge(dates, ap_charge, ap_charge_err, "")
    plot_ap_charge(dates, he_ap_charge, he_ap_charge_err, "he")

    plot_aapc_vs_charge(dates, ap_ratio, ap_ratio_err, "")
    plot_aapc_vs_charge(dates, he_ap_ratio, he_ap_ratio_err, "he")


if __name__ == "__main__":
    main()
