import sys

sys.path.insert(1, '../..')

# import ROOT and bash commands
import ROOT
import tqdm

# import python plotting and numpy modules
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from array import array
from scipy.signal import find_peaks

# import stats module
from scipy.optimize import curve_fit

# import custom made classes
from functions.other_functions import pmt_parse_arguments, fit, chi2, process_date, linear, gaus
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


def create_file(file_name: str):
    file = open(file_name, 'w')
    file.close()


def write_to_file(file_name: str, line):
    file = open(file_name, 'a')

    file.write(line+'\n')

    file.close()


def get_error_a_divide_b(da, a, db, b):
    c = a/b
    dc = c * np.sqrt((da/a)**2 + (db/b)**2)
    return dc


def get_resolution(mu: float, mu_err: float, sig: float, sig_err: float):
    res = 0
    res_err = 0
    if mu == 0:
        pass
    else:
        res = sig/mu
        res_err = get_error_a_divide_b(sig_err, sig, mu_err, mu)
    return res*100, res_err*100


def read_file(date: str, voltage: int, root_file_name: str, pmt_array: PMT_Array, output_file_location: str):

    file = ROOT.TFile(root_file_name, "READ")
    file.cd()

    apulse_info = [[] for i in range(pmt_array.get_pmt_total_number())]

    for i_om in range(pmt_array.get_pmt_total_number()):
        apulse_num_hist = file.Get(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() +
                                   "_apulse_num_" + str(voltage) + "V")
        apulse_time_hist = file.Get(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() +
                                    "_apulse_times_" + str(voltage) + "V")
        apulse_amplitude_hist = file.Get(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() +
                                         "_apulse_amplitudes_" + str(voltage) + "V")

        he_apulse_num_hist = file.Get(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() +
                                      "_he_apulse_num_" + str(voltage) + "V")
        he_apulse_amplitude_hist = file.Get(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() +
                                            "_he_apulse_amplitudes_" + str(voltage) + "V")
        ap_charge_hist = file.Get(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() +
                                  "_ap_charge_spectrum_" + str(voltage) + "V")
        he_ap_charge_hist = file.Get(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() +
                                     "_he_ap_charge_spectrum_" + str(voltage) + "V")
        ap_charge_charge_hist = file.Get(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() +
                                         "_ap_charge_charge_spectrum_" + str(voltage) + "V")
        he_ap_charge_charge_hist = file.Get(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() +
                                            "_he_ap_charge_charge_spectrum_" + str(voltage) + "V")

        try:
            apulse_num_hist.GetEntries()
            apulse_time_hist.GetEntries()
            apulse_amplitude_hist.GetEntries()
            he_apulse_num_hist.GetEntries()
            he_apulse_amplitude_hist.GetEntries()
            ap_charge_hist.GetEntries()
            he_ap_charge_hist.GetEntries()
        except:
            continue

        apulse_rate = 0
        for i_bin in range(2, apulse_num_hist.GetNbinsX()):
            apulse_rate += apulse_num_hist.GetBinContent(i_bin)

        par = (apulse_rate/apulse_num_hist.GetEntries()) * 100
        par_err = (np.sqrt(1/apulse_rate + 1/apulse_num_hist.GetEntries())) * apulse_rate

        he_apulse_rate = 0
        for i_bin in range(2, he_apulse_num_hist.GetNbinsX()):
            he_apulse_rate += he_apulse_num_hist.GetBinContent(i_bin)

        par_he = (he_apulse_rate/apulse_num_hist.GetEntries()) * 100
        par_he_err = (np.sqrt(1/he_apulse_rate + 1/apulse_num_hist.GetEntries())) * he_apulse_rate

        aan = apulse_num_hist.GetMean()
        aan_err = apulse_num_hist.GetMeanError()
        aan_he = he_apulse_num_hist.GetMean()
        aan_he_err = he_apulse_num_hist.GetMeanError()
        ap_charge = ap_charge_hist.GetMean()
        ap_charge_err = ap_charge_hist.GetMeanError()
        he_ap_charge = he_ap_charge_hist.GetMean()
        he_ap_charge_err = he_ap_charge_hist.GetMeanError()
        ratio = ap_charge_charge_hist.GetMean()
        he_ratio = he_ap_charge_charge_hist.GetMean()

        pars = {
            "par": par,
            "par_err": par_err,
            "par_he": par_he,
            "par_he_err": par_he_err,
            "aan": aan,
            "aan_err": aan_err,
            "aan_he": aan_he,
            "aan_he_err": aan_he_err,
            "ap_charge": ap_charge,
            "ap_charge_err": ap_charge_err,
            "he_ap_charge": he_ap_charge,
            "he_ap_charge_err": he_ap_charge_err,
            "ratio": ratio,
            "he_ratio": he_ratio
        }
        apulse_info[i_om].append(pars)

    file.Close()
    return apulse_info


def plot_par(date, par, output_directory: str, pmt_object: PMT_Object, name: str):
    date = process_date(date)
    try:
        start = np.where(date == 0)[0][0]
    except:
        start = np.where(date == 1)[0][0]
    mid = np.where(date == 98)[0][0]

    popt, pcov = curve_fit(f=model_day, xdata=date[mid + 1:], ydata=np.array(par[mid + 1:]),
                           p0=[0.1, 1, 19], bounds=[[0, 0, 0], [1e5, 1e10, 100]])

    plt.figure(num=None, figsize=(9, 5), dpi=80, facecolor='w', edgecolor='k')
    plt.plot(date[:start + 1], np.array(par[:start + 1]), "g.", label="Atmospheric He")
    plt.plot(date[start + 1:mid + 1], np.array(par[start + 1:mid + 1]), "b.", label="1% He")
    plt.plot(date[mid + 1:], np.array(par[mid + 1:]), "r.", label="10% He")
    plt.plot(date[mid + 1:], model_day(date[mid + 1:], *popt), 'k-', label='model')
    plt.axvline(date[start], 0, 100, ls='--', color='k')
    plt.axvline(date[mid], 0, 100, ls='--', color='k')
    plt.xlabel("exposure days relative to 06/11/2019")
    plt.ylabel("Afterpulse rate /%")
    plt.title(pmt_object.get_pmt_id() + " PAR vs exposure time")
    plt.grid()
    plt.ylim(0, 100)
    plt.legend(loc='upper left')
    plt.savefig(output_directory + "/summary_plots/" +
                pmt_object.get_pmt_id() + "_par_vs_time" + name + ".pdf")

    print('Param p{}: {:.2e} ± {:.1e}'.format(0, popt[0], np.sqrt(pcov[0, 0])))
    print('Param p{}: {:.2e} ± {:.1e}'.format(1, popt[1], np.sqrt(pcov[1, 1])))
    print('Param p{}: {:.2e} ± {:.1e}'.format(2, popt[2], np.sqrt(pcov[2, 2])))
    # print('Chi2 is:', chi_2)

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


def plot_aan(dates, aans, output_directory: str, name: str):
    date_0 = process_date(dates[0])
    date_1 = process_date(dates[1])
    try:
        start = np.where(date_0 == 0)[0][0]
    except:
        start = np.where(date_0 == 1)[0][0]
    mid = np.where(date_0 == 98)[0][0]

    '''popt, pcov = curve_fit(f=model_day, xdata=date[mid + 1:], ydata=np.array(aan[mid + 1:]),
                           p0=[0.01, 0.01, 0.8], bounds=[[0, 0, 0], [1e5, 1e10, 100]])'''
    aan_0 = np.array(aans[0])
    aan_1 = np.array(aans[1])

    x = np.array(date_0[mid + 1:]) - date_0[mid + 1:][0]
    y = aan_0[mid + 1:]
    yerr = aan_0[mid + 1:]*0.01
    pars, errs, chi = root_fit(x, y, yerr)

    fig1 = plt.figure(num=None, figsize=(9, 5), dpi=80, facecolor='w', edgecolor='k')
    frame1 = fig1.add_axes((.1, .3, .8, .6))
    frame1.set_xticklabels([])
    plt.errorbar(date_0[:start + 1], aan_0[:start + 1], zorder=0,
                 yerr=aan_0[:start + 1]*0.01, fmt="g.", label="Atmospheric He")
    plt.errorbar(date_0[start + 1:mid + 1], aan_0[start + 1:mid + 1], zorder=0,
                 yerr=aan_0[start + 1:mid + 1]*0.01, fmt="b.", label="1% He")
    plt.errorbar(date_0[mid + 1:], aan_0[mid + 1:], zorder=0,
                 yerr=aan_0[mid + 1:]*0.01, fmt="r.", label="10% He")
    plt.errorbar(date_1, aan_1, zorder=0,
                 yerr=aan_1*0.01, fmt="C1.", label="Control")
    plt.plot(date_0[mid + 1:], Model().func(x, pars), 'k-', label='model', zorder=10,)
    plt.axvline(date_0[start], 0, 100, ls='--', color='k')
    plt.axvline(date_0[mid], 0, 100, ls='--', color='k')

    handles, labels = plt.gca().get_legend_handles_labels()
    patch = matplotlib.patches.Patch(color='white', label=r'$P_0 =$ {:.4e} ± {:.0e}'.format(pars[0], errs[0]))
    patch_1 = matplotlib.patches.Patch(color='white', label=r'$P_1 =$ {:.4e} ± {:.0e}'.format(pars[1], errs[1]))
    patch_2 = matplotlib.patches.Patch(color='white', label=r'$L =$ {:.0f} ± {:.0f}'.format(pars[2]/(3600*24), errs[2]/(3600*24)))
    patch_3 = matplotlib.patches.Patch(color='white', label=r'$\chi^2_R =$ {:.2f}'.format(chi))
    handles.extend([patch, patch_1, patch_2, patch_3])

    plt.ylabel("Average afterpulse number")
    plt.title("AAN vs exposure time")
    plt.grid()
    plt.legend(handles=handles, loc='upper left')
    plt.xlim(-30, 420)

    frame2 = fig1.add_axes((.1, .1, .8, .2))
    plt.xlabel("exposure days relative to 06/11/2019")
    plt.axhline(0, ls='--', color='black')
    plt.ylabel("(model-data)/model")
    plt.grid()
    plt.xlim(-30, 420)
    plt.errorbar(date_0[mid + 1:],
                 (Model().func(x, pars) - aan_0[mid + 1:])/Model().func(x, pars),
                 yerr=aan_0[mid + 1:]*0.01/Model().func(date_0[mid + 1:], pars), fmt="k.")

    plt.savefig(output_directory + "/summary_plots/" + "aan_vs_time" + name + ".pdf")
    plt.close()


def plot_par_ratio(dates: list, par: list, output_directory: str, name: str):
    # Plot ratio
    x_date = []
    ratio = []
    ratio_err = []
    for i in range(len(dates[0])):
        for j in range(len(dates[1])):
            if dates[0][i] == dates[1][j]:
                x_date.append(dates[0][i])
                if par[1][j] == 0:
                    pass
                else:
                    ratio.append(par[0][i] / par[1][j])
                break

    x_date = process_date(x_date)

    plt.plot(x_date, ratio, "k.")
    plt.axvline(98, color="r", ls="--")
    plt.axvline(0, color="b", ls="--")
    plt.xlabel("exposure days relative to 06/11/2019")
    plt.ylabel("Ratio PAR Ch0/Ch1")
    plt.title("Ratio of PAR vs time")
    plt.grid()
    # plt.xlim(np.amin(np.array(x_date)), np.amax(np.array(x_date)))
    # plt.ylim(0, 2)
    plt.savefig(output_directory + "/summary_plots/par_ratio_vs_time" + name + ".pdf")
    plt.close()


def plot_aan_ratio(dates: list, aan: list, output_directory: str, name: str):
    # Plot ratio
    x_date = []
    ratio = []
    ratio_err = []
    for i in range(len(dates[0])):
        for j in range(len(dates[1])):
            if dates[0][i] == dates[1][j]:
                x_date.append(dates[0][i])
                if aan[1][j] == 0:
                    pass
                else:
                    ratio.append(aan[0][i] / aan[1][j])
                break

    x_date = process_date(x_date)

    plt.plot(x_date, ratio, "k.")
    plt.axvline(98, color="r", ls="--")
    plt.axvline(0, color="b", ls="--")
    plt.xlabel("exposure days relative to 06/11/2019")
    plt.ylabel("Ratio AAN Ch0/Ch1")
    plt.title("Ratio of AAN vs time")
    plt.grid()
    # plt.xlim(np.amin(np.array(x_date)), np.amax(np.array(x_date)))
    # plt.ylim(0, 2)
    plt.savefig(output_directory + "/summary_plots/aan_ratio_vs_time" + name + ".pdf")
    plt.close()


def plot_ap_charge(dates, ap_charges, output_directory: str, name: str):
    # print(len(ap_charge), ap_charge)
    date_0 = process_date(dates[0])
    date_1 = process_date(dates[1])
    try:
        start = np.where(date_0 == 0)[0][0]
    except:
        start = np.where(date_0 == 1)[0][0]
    mid = np.where(date_0 == 98)[0][0]

    ap_charge_0 = np.array(ap_charges[0])
    ap_charge_1 = np.array(ap_charges[1])

    x = np.array(date_0[mid + 1:]) - date_0[mid + 1:][0]
    y = np.array(ap_charge_0[mid + 1:])
    yerr = np.array(ap_charge_0[mid + 1:]) * 0.01
    pars, errs, chi = root_fit(x, y, yerr)

    fig1 = plt.figure(num=None, figsize=(9, 5), dpi=80, facecolor='w', edgecolor='k')
    frame1 = fig1.add_axes((.1, .3, .8, .6))
    frame1.set_xticklabels([])
    plt.errorbar(date_0[:start + 1], ap_charge_0[:start + 1], zorder=0,
                 yerr=ap_charge_0[:start + 1]*0.01, fmt="g.", label="Atmospheric He")
    plt.errorbar(date_0[start + 1:mid + 1], ap_charge_0[start + 1:mid + 1], zorder=0,
                 yerr=ap_charge_0[start + 1:mid + 1]*0.01, fmt="b.", label="1% He")
    plt.errorbar(date_0[mid + 1:], ap_charge_0[mid + 1:], zorder=0,
                 yerr=ap_charge_0[mid + 1:]*0.01, fmt="r.", label="10% He")
    plt.errorbar(date_1, ap_charge_1, zorder=0,
                 yerr=ap_charge_1*0.01, fmt="C1.", label="Control")
    plt.plot(date_0[mid + 1:], Model().func(x, pars), 'k-', label='model', zorder=10)
    plt.axvline(date_0[start], 0, 100, ls='--', color='k')
    plt.axvline(date_0[mid], 0, 100, ls='--', color='k')

    handles, labels = plt.gca().get_legend_handles_labels()
    patch = matplotlib.patches.Patch(color='white', label=r'$P_0 =$ {:.4e} ± {:.0e}'.format(pars[0], errs[0]))
    patch_1 = matplotlib.patches.Patch(color='white', label=r'$P_1 =$ {:.4e} ± {:.0e}'.format(pars[1], errs[1]))
    patch_2 = matplotlib.patches.Patch(color='white', label=r'$L =$ {:.0f} ± {:.0f}'.format(pars[2] / (3600 * 24),
                                                                                            errs[2] / (3600 * 24)))
    patch_3 = matplotlib.patches.Patch(color='white', label=r'$\chi^2_R =$ {:.2f}'.format(chi))
    handles.extend([patch, patch_1, patch_2, patch_3])

    plt.ylabel("Average Charge /pC")
    plt.title("Average After-pulse Region charge")
    plt.grid()
    plt.xlim(-30, 420)
    plt.legend(handles=handles, loc='upper left')

    frame2 = fig1.add_axes((.1, .1, .8, .2))
    plt.xlabel("exposure days relative to 06/11/2019")
    plt.axhline(0, ls='--', color='black')
    plt.ylabel("(model-data)/model")
    plt.grid()
    plt.xlim(-30, 420)
    plt.errorbar(date_0[mid + 1:],
                 (Model().func(x, pars) - ap_charge_0[mid + 1:]) / Model().func(x, pars),
                 yerr=ap_charge_0[mid + 1:] * 0.01 / Model().func(date_0[mid + 1:], pars), fmt="k.")

    plt.savefig(output_directory + "/summary_plots/" + "ap_charge_vs_time" + name + ".pdf")
    plt.close()


def plot_aapc_vs_charge(dates, ratios, output_directory: str, name: str):
    date_0 = process_date(dates[0])
    date_1 = process_date(dates[1])
    try:
        start = np.where(date_0 == 0)[0][0]
    except:
        start = np.where(date_0 == 1)[0][0]
    mid = np.where(date_0 == 98)[0][0]

    ratio_0 = np.array(ratios[0])
    ratio_1 = np.array(ratios[1])

    x = np.array(date_0[mid + 1:]) - date_0[mid + 1:][0]
    y = np.array(ratio_0[mid + 1:])
    yerr = np.array(ratio_0[mid + 1:]) * 0.01
    pars, errs, chi = root_fit(x, y, yerr)

    fig1 = plt.figure(num=None, figsize=(9, 5), dpi=80, facecolor='w', edgecolor='k')
    frame1 = fig1.add_axes((.1, .3, .8, .6))
    frame1.set_xticklabels([])
    plt.errorbar(date_0[:start + 1], ratio_0[:start + 1], zorder=0,
                 yerr=ratio_0[:start + 1]*0.01, fmt="g.", label="Atmospheric He")
    plt.errorbar(date_0[start + 1:mid + 1], ratio_0[start + 1:mid + 1], zorder=0,
                 yerr=ratio_0[start + 1:mid + 1]*0.01, fmt="b.", label="1% He")
    plt.errorbar(date_0[mid + 1:], ratio_0[mid + 1:], zorder=0,
                 yerr=ratio_0[mid + 1:]*0.01, fmt="r.", label="10% He")
    plt.errorbar(date_1, ratio_1, zorder=0,
                 yerr=ratio_1*0.01, fmt="C1.", label="Control")
    plt.plot(date_0[mid + 1:], Model().func(x, pars), 'k-', label='model', zorder=10)
    plt.axvline(date_0[start], 0, 100, ls='--', color='k')
    plt.axvline(date_0[mid], 0, 100, ls='--', color='k')

    handles, labels = plt.gca().get_legend_handles_labels()
    patch = matplotlib.patches.Patch(color='white', label=r'$P_0 =$ {:.4e} ± {:.0e}'.format(pars[0], errs[0]))
    patch_1 = matplotlib.patches.Patch(color='white', label=r'$P_1 =$ {:.4e} ± {:.0e}'.format(pars[1], errs[1]))
    patch_2 = matplotlib.patches.Patch(color='white', label=r'$L =$ {:.0f} ± {:.0f}'.format(pars[2] / (3600 * 24),
                                                                                            errs[2] / (3600 * 24)))
    patch_3 = matplotlib.patches.Patch(color='white', label=r'$\chi^2_R =$ {:.2f}'.format(chi))
    handles.extend([patch, patch_1, patch_2, patch_3])

    plt.ylabel("charge ratio")
    plt.title("Pulse charge apulse charge ratio")
    plt.grid()
    plt.legend(handles=handles, loc='upper left')
    plt.xlim(-30, 420)

    frame2 = fig1.add_axes((.1, .1, .8, .2))
    plt.xlabel("exposure days relative to 06/11/2019")
    plt.axhline(0, ls='--', color='black')
    plt.ylabel("(model-data)/model")
    plt.grid()
    plt.xlim(-30, 420)
    plt.errorbar(date_0[mid + 1:],
                 (Model().func(x, pars) - ratio_0[mid + 1:]) / Model().func(x, pars),
                 yerr=ratio_0[mid + 1:] * 0.01 / Model().func(date_0[mid + 1:], pars), fmt="k.")

    plt.savefig(output_directory + "/summary_plots/" + "ap_charge_charge_vs_time" + name + ".pdf")
    plt.close()


def main():
    # Handle the input arguments:
    ##############################
    args = pmt_parse_arguments()
    input_directory = args.i
    run_id = args.c
    # config_file_name = args.c
    output_directory = args.o
    ##############################

    filenames_txt = input_directory + "/filenames.txt"

    try:
        print(">>> Reading data from file: {}".format(filenames_txt))
        date_file = open(filenames_txt, 'r')
    except FileNotFoundError as fnf_error:
        print(fnf_error)
        raise Exception("Error opening data file {}".format(filenames_txt))

    filenames = np.loadtxt(filenames_txt, delimiter=',', dtype={
        'names': ['filename'],
        'formats': ['S100']}, unpack=True)

    topology = [2, 1]
    pmt_array = PMT_Array(topology, run_id+"summary")
    pmt_array.set_pmt_id("GAO607", 0)
    pmt_array.set_pmt_id("GAO612", 1)

    # Set up the containers for the summary
    par = [[] for i in range(pmt_array.get_pmt_total_number())]
    par_err = [[] for i in range(pmt_array.get_pmt_total_number())]
    par_he = [[] for i in range(pmt_array.get_pmt_total_number())]
    par_he_err = [[] for i in range(pmt_array.get_pmt_total_number())]
    aan = [[] for i in range(pmt_array.get_pmt_total_number())]
    aan_err = [[] for i in range(pmt_array.get_pmt_total_number())]
    aan_he = [[] for i in range(pmt_array.get_pmt_total_number())]
    aan_he_err = [[] for i in range(pmt_array.get_pmt_total_number())]
    ap_charge = [[] for i in range(pmt_array.get_pmt_total_number())]
    ap_charge_err = [[] for i in range(pmt_array.get_pmt_total_number())]
    he_ap_charge = [[] for i in range(pmt_array.get_pmt_total_number())]
    he_ap_charge_err = [[] for i in range(pmt_array.get_pmt_total_number())]
    ratios = [[] for i in range(pmt_array.get_pmt_total_number())]
    he_ratios = [[] for i in range(pmt_array.get_pmt_total_number())]
    dates = [[] for i in range(pmt_array.get_pmt_total_number())]

    for i_file in tqdm.tqdm(range(filenames.size)):
        file = filenames[i_file][0].decode("utf-8")

        date = file.split("_")[0]
        voltage = int(file.split("_")[1].split("A")[1])

        if voltage == 1400:
            pass
        else:
            continue

        apulse_info = read_file(date, voltage, input_directory + "/" + file, pmt_array, output_directory)

        for i_om in range(pmt_array.get_pmt_total_number()):

            if len(apulse_info[i_om]) > 0:
                pass
            else:
                continue

            i_par = apulse_info[i_om][0]["par"]
            i_par_he = apulse_info[i_om][0]["par_he"]
            i_par_err = apulse_info[i_om][0]["par_err"]
            i_par_he_err = apulse_info[i_om][0]["par_he_err"]

            i_aan = apulse_info[i_om][0]["aan"]
            i_aan_he = apulse_info[i_om][0]["aan_he"]
            i_aan_err = apulse_info[i_om][0]["aan_err"]
            i_aan_he_err = apulse_info[i_om][0]["aan_he_err"]
            i_ap_charge = apulse_info[i_om][0]["ap_charge"]
            i_ap_charge_err = apulse_info[i_om][0]["ap_charge_err"]
            i_he_ap_charge = apulse_info[i_om][0]["he_ap_charge"]
            i_he_ap_charge_err = apulse_info[i_om][0]["he_ap_charge_err"]
            i_ratio = apulse_info[i_om][0]["ratio"]
            i_he_ratio = apulse_info[i_om][0]["he_ratio"]

            par[i_om].append(i_par)
            par_err[i_om].append(i_par_err/10)
            par_he[i_om].append(i_par_he)
            par_he_err[i_om].append(i_par_he_err/10)
            aan[i_om].append(i_aan)
            aan_err[i_om].append(i_aan_err)
            aan_he[i_om].append(i_aan_he)
            aan_he_err[i_om].append(i_aan_he_err)
            ap_charge[i_om].append(i_ap_charge)
            ap_charge_err[i_om].append(i_ap_charge_err)
            he_ap_charge[i_om].append(i_he_ap_charge)
            he_ap_charge_err[i_om].append(i_he_ap_charge_err)
            ratios[i_om].append(i_ratio)
            he_ratios[i_om].append(i_he_ratio)

            dates[i_om].append(int(date))

    plot_ap_charge(dates, ap_charge, output_directory, "_" + run_id)
    plot_ap_charge(dates, he_ap_charge, output_directory, "_he_" + run_id)
    plot_aan(dates, aan, output_directory, "_" + run_id)
    plot_aan(dates, aan_he, output_directory, "_he_" + run_id)
    plot_aapc_vs_charge(dates, ratios, output_directory, "_" + run_id)
    plot_aapc_vs_charge(dates, he_ratios, output_directory, "_he_" + run_id)

    '''plot_par_ratio(dates, par, output_directory, "_" + run_id)
    plot_par_ratio(dates, par_he, output_directory, "_he_" + run_id)
    plot_aan_ratio(dates, aan, output_directory, "_" + run_id)
    plot_aan_ratio(dates, aan_he, output_directory, "_he_" + run_id)'''


if __name__ == '__main__':
    main()
