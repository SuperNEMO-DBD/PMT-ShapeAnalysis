import sys

sys.path.insert(1, '../..')

# import ROOT and bash commands
import ROOT
import tqdm

# import python plotting and numpy modules
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

# import stats module
from scipy.optimize import curve_fit

# import custom made classes
from functions.other_functions import pmt_parse_arguments, fit, chi2, process_date, linear, gaus
from src.PMT_Classes import *


def model_day(x, p0, p1, p2):
    y = []
    for i in range(len(x)):
        t = x[i]*3600*24
        temp = 0
        for n in range(1,11):
            temp += ((-1)**n/n**2)*(1 - np.exp(-(n**2)*(np.pi**2)*t/p1))
        f2 = (2/np.pi**2)*p1*temp
        y.append(p0*(t + f2) + p2)
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

        try:
            apulse_num_hist.GetEntries()
            apulse_time_hist.GetEntries()
            apulse_amplitude_hist.GetEntries()
            he_apulse_num_hist.GetEntries()
            he_apulse_amplitude_hist.GetEntries()
            ap_charge_hist.GetEntries()
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
            "ap_charge_err": ap_charge_err
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


def plot_aan(date, aan, output_directory: str, pmt_object: PMT_Object, name: str):
    date = process_date(date)
    try:
        start = np.where(date == 0)[0][0]
    except:
        start = np.where(date == 1)[0][0]
    mid = np.where(date == 98)[0][0]

    popt, pcov = curve_fit(f=model_day, xdata=date[mid + 1:], ydata=np.array(aan[mid + 1:]),
                           p0=[0.01, 0.01, 0.8], bounds=[[0, 0, 0], [1e5, 1e10, 100]])

    plt.figure(num=None, figsize=(9, 5), dpi=80, facecolor='w', edgecolor='k')
    plt.plot(date[:start + 1], np.array(aan[:start + 1]), "g.", label="Atmospheric He")
    plt.plot(date[start + 1:mid + 1], np.array(aan[start + 1:mid + 1]), "b.", label="1% He")
    plt.plot(date[mid + 1:], np.array(aan[mid + 1:]), "r.", label="10% He")
    plt.plot(date[mid + 1:], model_day(date[mid + 1:], *popt), 'k-', label='model')
    plt.axvline(date[start], 0, 100, ls='--', color='k')
    plt.axvline(date[mid], 0, 100, ls='--', color='k')
    plt.xlabel("exposure days relative to 06/11/2019")
    plt.ylabel("Average afterpulse number")
    plt.title(pmt_object.get_pmt_id() + " AAN vs exposure time")
    plt.grid()
    # plt.ylim(0, 1.5)
    plt.legend(loc='upper left')
    plt.savefig(output_directory + "/summary_plots/" +
                pmt_object.get_pmt_id() + "_aan_vs_time" + name + ".pdf")
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


def plot_ap_charge(date, ap_charge, output_directory: str, pmt_object: PMT_Object, name: str):
    # print(len(ap_charge), ap_charge)
    date = process_date(date)
    try:
        start = np.where(date == 0)[0][0]
    except:
        start = np.where(date == 1)[0][0]
    mid = np.where(date == 98)[0][0]

    plt.figure(num=None, figsize=(9, 5), dpi=80, facecolor='w', edgecolor='k')
    plt.plot(date[:start + 1], np.array(ap_charge[:start + 1]), "g.", label="Atmospheric He")
    plt.plot(date[start + 1:mid + 1], np.array(ap_charge[start + 1:mid + 1]), "b.", label="1% He")
    plt.plot(date[mid + 1:], np.array(ap_charge[mid + 1:]), "r.", label="10% He")
    plt.axvline(date[start], 0, 100, ls='--', color='k')
    plt.axvline(date[mid], 0, 100, ls='--', color='k')
    plt.xlabel("exposure days relative to 06/11/2019")
    plt.ylabel("Average Charge /pC")
    plt.title(pmt_object.get_pmt_id() + " Average After-pulse Region charge")
    plt.grid()
    # plt.ylim(0, 1.5)
    plt.legend(loc='upper left')
    plt.savefig(output_directory + "/summary_plots/" +
                pmt_object.get_pmt_id() + "_ap_charge_vs_time" + name + ".pdf")
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

            dates[i_om].append(int(date))

    for i_om in range(pmt_array.get_pmt_total_number()):

        plot_par(dates[i_om], par[i_om], output_directory, pmt_array.get_pmt_object_number(i_om), "_" + run_id)
        plot_par(dates[i_om], par_he[i_om], output_directory, pmt_array.get_pmt_object_number(i_om), "_he_" + run_id)
        plot_aan(dates[i_om], aan[i_om], output_directory, pmt_array.get_pmt_object_number(i_om), "_" + run_id)
        plot_aan(dates[i_om], aan_he[i_om], output_directory, pmt_array.get_pmt_object_number(i_om), "_he_" + run_id)
        plot_ap_charge(dates[i_om], ap_charge[i_om], output_directory, pmt_array.get_pmt_object_number(i_om), "_" + run_id)

    plot_par_ratio(dates, par, output_directory, "_" + run_id)
    plot_par_ratio(dates, par_he, output_directory, "_he_" + run_id)
    plot_aan_ratio(dates, aan, output_directory, "_" + run_id)
    plot_aan_ratio(dates, aan_he, output_directory, "_he_" + run_id)


if __name__ == '__main__':
    main()
