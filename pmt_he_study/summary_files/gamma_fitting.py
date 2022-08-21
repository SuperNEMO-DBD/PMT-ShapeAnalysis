import sys

from xml.dom import minidom
from time import time
import json

import matplotlib.pyplot as plt
import numpy as np
import tqdm

sys.path.insert(1, '../..')
from pmt_he_study.models import *
from ReadRED import sndisplay as sn


def double_compton_spec(x, pars):
    return pars[0]*((1-1/(1+np.exp(-(x-pars[1])/pars[2]))) + gaussian_noh(x, pars[1], pars[4], pars[3])) +\
           pars[5]*((1-1/(1+np.exp(-(x-pars[6])/pars[7]))) + gaussian_noh(x, pars[6], pars[9], pars[8]))


def fill_spec(x, g_results, c_results):
    g_pars = g_results["pars"]
    c_pars = c_results["pars"]
    val = c_pars[0]*(g_pars[0]*((1-1/(1+np.exp(-(x-g_pars[1])/g_pars[2]))) + gaussian_noh(x, g_pars[1], g_pars[4], g_pars[3])) +
                     g_pars[5]*((1-1/(1+np.exp(-(x-g_pars[6])/g_pars[7]))) + gaussian_noh(x, g_pars[6], g_pars[9], g_pars[8]))) +\
          c_pars[1]*(gaussian_noh(x, c_pars[2], c_pars[3], 1.54) + gaussian_noh(x, c_pars[2]*(1 + 72.144/481.694), c_pars[3]*1.072, 0.44) + gaussian_noh(x, c_pars[2]*(1 + 84.154/481.694), c_pars[3]*1.084, 0.11)) +\
          c_pars[4]*(gaussian_noh(x, c_pars[5], c_pars[6], 7.08) + gaussian_noh(x, c_pars[5]*(1 + 72.144/975.651), c_pars[6]*1.036, 1.84) + gaussian_noh(x, c_pars[5]*(1 + 84.154/975.651), c_pars[6]*1.042, 0.44))
    return val


def ce_and_comp_spec(x, g_results, c_results):
    g_pars = g_results["pars"]
    c_pars = c_results["pars"]
    val = 0.752515*(g_pars[0]*((1 - 1 / (1 + np.exp(-(x - g_pars[1]) / g_pars[2]))) + gaussian_noh(x, g_pars[1], g_pars[4], g_pars[3])) + g_pars[5] * ((1 - 1 / (1 + np.exp(-(x - g_pars[6]) / g_pars[7]))) + gaussian_noh(x, g_pars[6], g_pars[9], g_pars[8]))) + \
          c_pars[0]*(gaussian_noh(x, c_pars[1], c_pars[2], 7.08) + gaussian_noh(x, c_pars[1] * (1 + 72.144 / 975.651),c_pars[2] * 1.036, 1.84) + gaussian_noh(x, c_pars[1] * (1 + 84.154 / 975.651), c_pars[2] * 1.042, 0.44))
    return val


def compton_spec(x, A, loc, width, B, sig):
    return A * ((1 - 1 / (1 + np.exp(-(x - loc) / width))) + gaussian_noh(x, loc, sig, B))


def triple_gaus(x, A, mu, sig):
    y = A * (gaussian_noh(x, mu, sig, 7.08) + gaussian_noh(x, mu * (1 + 72.144 / 975.651),sig * 1.036, 1.84) + gaussian_noh(x,mu * (1 + 84.154 / 975.651),sig * 1.042,0.44))
    return y


def root_fit_gamma(data, plot=False):
    """
    This function plots the spectrum from the gamma interactions. This includes 2 Compton edges and
    2 photo peaks (gaussian)
    :param data: data from json file loading
    :param plot: boolean for plotting
    :return: ROOT TMinuit fit results dict
    """
    ######################################################
    #   - The function for a compton edge: 1 - 1/(1 + exp(-(x - A)/B))
    #            where A is the compton edge energy and B is the slope
    #   - the equation for A is from this formula
    #       - E_c = E(1 - 1/(1 + 2E/m)) (m in keV = 511keV)
    #   - E1 = 569.698 keV and E2 = 1063.656 keV
    #
    #
    ######################################################
    hist = ROOT.TH1D("gamma_fit", "gamma_fit", 240, 0, 60)
    for i in range(len(data['charge'])):
        hist.Fill(data['charge'][i])

    lower, higher = 12, 45
    x = np.arange(lower, higher, 0.25)
    fit = ROOT.TF1("fit",
                   "[0]*((1-1/(1+TMath::Exp(-(x-[1])/[2])))"
                   "+"
                   "[3]*TMath::Gaus(x,[1],[4]))"
                   "+"
                   "[5]*((1-1/(1+TMath::Exp(-(x-[6])/[7])))"
                   "+"
                   "[8]*TMath::Gaus(x,[6],[9]))",
                   lower, higher)

    names = ["A_0", "loc_0", "width_0", "B_0", "sig_0",
             "A_1", "loc_1", "width_1", "B_1", "sig_1"]
    guess = [1023, 15, 1.5, 620, 1.25,
             500, 40, 1.5, 500, 1.25]
    bounds = [
        [0, lower, 0, 0, 1,  0, 20, 0, 0, 1],
        [10000, 20, 100, 10000, 100, 10000, higher, 100, 10000, 100]
    ]
    fit.SetParNames(*names)
    fit.SetParameters(*guess)
    for i in range(len(names)):
        fit.SetParLimits(i, bounds[0][i], bounds[1][i])
    hist.Fit("fit", "", "", lower, higher)

    pars = []
    fit_result = {
        "chi2": fit.GetChisquare(),
        "ndof": fit.GetNDF()
    }
    for i in range(len(names)):
        fit_result[names[i]] = (fit.GetParameter(i), fit.GetParError(i))
        pars.append(fit.GetParameter(i))
    fit_result["pars"] = pars

    if plot:
        freq, bin_edges = np.histogram(data['charge'], range=(0, 60), bins=240)
        width = bin_edges[2] - bin_edges[1]
        bin_centres = bin_edges[:-1] + width / 2

        y_model = double_compton_spec(x, pars)
        y_data = freq[int(lower/width): int(higher/width)]
        y_data_err = np.sqrt(y_data)

        fig = plt.figure(figsize=(5, 5))
        frame1 = fig.add_axes((.15, .32, .8, .6))
        frame1.set_xticklabels([])

        plt.bar(x=bin_centres, height=freq, width=width, color="C1", alpha=0.5)
        plt.errorbar(x=bin_centres, y=freq, yerr=np.sqrt(freq), fmt='k.', label='data',
                     markersize=marker_size, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
        plt.plot(x, y_model, "r-",
                 label=r'model $\chi^2$/N$_{DoF}$:' + '{:.2f}/{}'.format(fit_result['chi2'], fit_result['ndof']))
        plt.plot(x, compton_spec(x, *pars[:5]), "C2--",
                 label=r'569.698keV')
        plt.plot(x, compton_spec(x, *pars[5:]), "C3--",
                 label=r'1063.656keV')
        plt.xlim(0, 60)
        plt.legend(loc='upper right')
        plt.ylabel("Counts")

        frame2 = fig.add_axes((.15, .1, .8, .2))
        plt.axhline(0, ls='--', color='k')
        plt.errorbar(x, (y_data - y_model)/y_model, yerr=y_data_err/y_model, fmt="k.",
                     markersize=marker_size, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
        plt.xlim(0, 60)
        plt.ylim(-0.25, 0.25)
        plt.xlabel("Charge /pC")
        plt.ylabel('(data-model)/data')

        plt.tight_layout()
        plt.savefig("/Users/williamquinn/Desktop/PMT_Project/gamma_fit.pdf")
        plt.close(fig)

    '''del hist
    del fit'''

    return fit_result


def root_fit_ce(data, gamma_results, plot=False):
    hist = ROOT.TH1D("full_fit", "full_fit", 240, 0, 60)
    for i in range(len(data['charge'])):
        hist.Fill(data['charge'][i])

    lower, higher = 17, 50
    x = np.arange(lower, higher, 0.25)
    g_pars = gamma_results["pars"]
    fit = ROOT.TF1("fit",
                   "[0]*({}*((1-1/(1+TMath::Exp(-(x-{})/{})))".format(g_pars[0], g_pars[1], g_pars[2]) +
                   "+" +
                   "{}*TMath::Gaus(x,{},{}))".format(g_pars[3], g_pars[1], g_pars[4]) +
                   "+" +
                   "{}*((1-1/(1+TMath::Exp(-(x-{})/{})))".format(g_pars[5], g_pars[6], g_pars[7]) +
                   "+" +
                   "{}*TMath::Gaus(x,{},{})))".format(g_pars[8], g_pars[6], g_pars[9]) +
                   " + [1]*(1.54*TMath::Gaus(x, [2], [3]) + 0.44*TMath::Gaus(x, [2]*(1 + 72.144/481.694), [3]*1.072) + 0.11*TMath::Gaus(x, [2]*(1 + 84.154/481.694), [3]*1.084))" +
                   " + [4]*(7.08*TMath::Gaus(x, [5], [6]) + 1.84*TMath::Gaus(x, [5]*(1 + 72.144/975.651), [6]*1.036) + 0.44*TMath::Gaus(x, [5]*(1 + 84.154/975.651), [6]*1.042))",
                   lower, higher)

    names = ["A", "A_0", "mu_0", "sig_0", "A_1", "mu_1", "sig_1"]
    guess = [100, 100, 20, 1, 100, 42, 1]
    bounds = [
        [0, 0, 10, 0, 0, 30, 0],
        [1000, 1000, 30, 10, 1000, 60, 10]
    ]
    fit.SetParNames(*names)
    fit.SetParameters(*guess)
    for i in range(len(names)):
        fit.SetParLimits(i, bounds[0][i], bounds[1][i])
    hist.Fit("fit", "", "", lower, higher)

    pars = []
    fit_result = {
        "chi2": fit.GetChisquare(),
        "ndof": fit.GetNDF()
    }
    for i in range(len(names)):
        fit_result[names[i]] = (fit.GetParameter(i), fit.GetParError(i))
        pars.append(fit.GetParameter(i))
    fit_result["pars"] = pars

    print(">>", fit_result["pars"], gamma_results["pars"])

    if plot:
        freq, bin_edges = np.histogram(data['charge'], range=(0, 60), bins=240)
        width = bin_edges[2] - bin_edges[1]
        bin_centres = bin_edges[:-1] + width / 2

        y_model = fill_spec(x, gamma_results, fit_result)
        y_data = freq[int(lower / width): int(higher / width)]
        y_data_err = np.sqrt(y_data)

        fig = plt.figure(figsize=(5, 5))
        frame1 = fig.add_axes((.15, .32, .8, .6))
        frame1.set_xticklabels([])

        plt.bar(x=bin_centres, height=freq, width=width, color="C3", alpha=0.5)
        plt.errorbar(x=bin_centres, y=freq, yerr=np.sqrt(freq), fmt='k.', label='data',
                     markersize=marker_size, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
        plt.plot(x, y_model, "r-",
                 label=r'model $\chi^2$/N$_{DoF}$:' + '{:.2f}/{}'.format(fit_result['chi2'], fit_result['ndof']))
        plt.xlim(0, 60)
        plt.legend(loc='upper left')
        plt.ylabel("Counts")

        frame2 = fig.add_axes((.15, .1, .8, .2))
        plt.axhline(0, ls='--', color='k')
        plt.errorbar(x, (y_data - y_model) / y_model, yerr=y_data_err / y_model, fmt="k.",
                     markersize=marker_size, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
        plt.xlim(0, 60)
        plt.ylim(-0.5, 0.5)
        plt.xlabel("Charge /pC")
        plt.ylabel('(data-model)/data')

        plt.tight_layout()
        plt.savefig("/Users/williamquinn/Desktop/PMT_Project/full_bi_spec_fit.pdf")
        plt.close(fig)

    del hist
    del fit

    return fit_result


def root_fit_both(data, gamma_results, plot=False):
    hist_0 = ROOT.TH1D("full_fit_bkgd", "full_fit_bkgd", 240, 0, 60)
    hist_1 = ROOT.TH1D("full_fit", "full_fit", 240, 0, 60)
    for i in range(len(data['charge'])):
        hist_0.Fill(data['charge'][i])
        hist_1.Fill(data['charge'][i])

    lower, higher = 41, 50
    x = np.arange(lower, higher, 0.25)
    g_pars = gamma_results["pars"]
    names = ["A", "mu", "sig"]
    guess = [100, 43, 1]
    bounds = [
        [0, lower, 0],
        [1000, higher, 5]
    ]

    fit_0 = ROOT.TF1("fit_0",
                   "0.752515*({}*((1-1/(1+TMath::Exp(-(x-{})/{})))".format(g_pars[0], g_pars[1], g_pars[2]) +
                   "+" +
                   "{}*TMath::Gaus(x,{},{}))".format(g_pars[3], g_pars[1], g_pars[4]) +
                   "+" +
                   "{}*((1-1/(1+TMath::Exp(-(x-{})/{})))".format(g_pars[5], g_pars[6], g_pars[7]) +
                   "+" +
                   "{}*TMath::Gaus(x,{},{})))".format(g_pars[8], g_pars[6], g_pars[9]) +
                   " + [0]*(7.08*TMath::Gaus(x, [1], [2]) + 1.84*TMath::Gaus(x, [1]*(1 + 72.144/975.651), [2]*1.036) + 0.44*TMath::Gaus(x, [1]*(1 + 84.154/975.651), [2]*1.042))",
                   lower, higher)
    fit_0.SetParNames(*names)
    fit_0.SetParameters(*guess)
    for i in range(len(names)):
        fit_0.SetParLimits(i, bounds[0][i], bounds[1][i])
    hist_0.Fit("fit_0", "", "", lower, higher)

    pars = []
    fit_result_0 = {
        "chi2": fit_0.GetChisquare(),
        "ndof": fit_0.GetNDF()
    }
    for i in range(len(names)):
        fit_result_0[names[i]] = (fit_0.GetParameter(i), fit_0.GetParError(i))
        pars.append(fit_0.GetParameter(i))
    fit_result_0["pars"] = pars

    fit_1 = ROOT.TF1("fit_1",
                     "[0]*(7.08*TMath::Gaus(x, [1], [2]) + 1.84*TMath::Gaus(x, [1]*(1 + 72.144/975.651), [2]*1.036) + 0.44*TMath::Gaus(x, [1]*(1 + 84.154/975.651), [2]*1.042))",
                     lower, higher)
    fit_1.SetParNames(*names)
    fit_1.SetParameters(*guess)
    for i in range(len(names)):
        fit_1.SetParLimits(i, bounds[0][i], bounds[1][i])
    hist_1.Fit("fit_1", "", "", lower, higher)

    pars = []
    fit_result_1 = {
        "chi2": fit_1.GetChisquare(),
        "ndof": fit_1.GetNDF()
    }
    for i in range(len(names)):
        fit_result_1[names[i]] = (fit_1.GetParameter(i), fit_1.GetParError(i))
        pars.append(fit_1.GetParameter(i))
    fit_result_1["pars"] = pars


    if plot:
        freq, bin_edges = np.histogram(data['charge'], range=(0, 60), bins=240)
        width = bin_edges[2] - bin_edges[1]
        bin_centres = bin_edges[:-1] + width / 2

        y_model_0 = ce_and_comp_spec(x, gamma_results, fit_result_0)
        y_model_1 = triple_gaus(x, fit_result_1["A"][0], fit_result_1["mu"][0], fit_result_1["sig"][0])

        res_0 = fit_result_0["sig"][0]/fit_result_0["mu"][0]
        res_0_err = res_0 * np.sqrt((fit_result_0["sig"][1]/fit_result_0["sig"][0])**2 + (fit_result_0["mu"][1]/fit_result_0["mu"][0])**2)
        res_1 = fit_result_1["sig"][0]/fit_result_1["mu"][0]
        res_1_err = res_1 * np.sqrt((fit_result_1["sig"][1]/fit_result_1["sig"][0])**2 + (fit_result_1["mu"][1]/fit_result_1["mu"][0])**2)

        res_0, res_1 = res_0*100, res_1*100
        res_0_err, res_1_err = res_0_err * 100, res_1_err * 100

        y_data = freq[int(lower / width): int(higher / width)]
        y_data_err = np.sqrt(y_data)

        fig = plt.figure(figsize=(5, 5))
        frame1 = fig.add_axes((.15, .32, .8, .6))
        frame1.set_xticklabels([])

        plt.bar(x=bin_centres, height=freq, width=width, color="C2", alpha=0.5)
        plt.errorbar(x=bin_centres, y=freq, yerr=np.sqrt(freq), fmt='k.', label='data',
                     markersize=marker_size, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
        plt.plot(x, y_model_0, "r-",
                 label=r'model 0 $\chi^2$/N$_{DoF}$:' +
                       '{:.2f}/{}'.format(fit_result_0['chi2'],
                                          fit_result_0['ndof']) + " R={:.1f}±{:.1f}%".format(res_0, res_0_err))
        plt.plot(x, y_model_1, "g--",
                 label=r'model 1 $\chi^2$/N$_{DoF}$:' +
                       '{:.2f}/{}'.format(fit_result_1['chi2'],
                                          fit_result_1['ndof']) + " R={:.1f}±{:.1f}%".format(res_1, res_1_err))
        plt.xlim(40, 52)
        plt.legend(loc='upper right')
        plt.ylabel("Counts")

        frame2 = fig.add_axes((.15, .1, .8, .2))
        plt.axhline(0, ls='--', color='k')
        plt.errorbar(x, (y_data - y_model_0) / y_model_0, yerr=y_data_err / y_model_0, fmt="r.",
                     markersize=marker_size, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
        plt.errorbar(x, (y_data - y_model_1) / y_model_1, yerr=y_data_err / y_model_1, fmt="g.",
                     markersize=marker_size, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
        plt.xlim(40, 52)
        plt.ylim(-0.5, 0.5)
        plt.xlabel("Charge /pC")
        plt.ylabel('(data-model)/data')

        plt.tight_layout()
        plt.savefig("/Users/williamquinn/Desktop/PMT_Project/bi_spec_comp.pdf")
        plt.close(fig)

    '''del hist_0
    del hist_1
    del fit_0
    del fit_1'''


def parse_data(input_data_file_name, name):
    print(">>> Parsing the data file...")
    processing_start = time()

    # parse an xml file by name
    xml_file = minidom.parse(input_data_file_name)
    events = xml_file.getElementsByTagName('event')
    parse_time = time() - processing_start

    print(">>> File is good. Parse time: %.3f s" % parse_time)
    print(">>> Number of Events: {}".format(len(events)))

    data = {"baselines": [], "amplitude": [], "charge": [], "peak": []}

    for event in tqdm.tqdm(range(len(events))):
        traces = events[event].getElementsByTagName('trace')
        for trace_index, trace in enumerate(traces):
            waveform = [int(i.strip()) for i in trace.firstChild.data.split(" ")[:-1]]
            baseline = float(np.average(waveform[:150]))
            peak = int(np.argmin(waveform))
            amplitude = float(waveform[peak]) - baseline
            charge = float(-1*(np.sum(waveform[peak - 10: peak+20]) - baseline*30)/50)
            data["baselines"].append(baseline)
            data["amplitude"].append(amplitude)
            data["peak"].append(peak)
            data["charge"].append(charge)
            # waveform_data_list.append(waveform)
    with open("/Users/williamquinn/Desktop/PMT_Project/" + name + "_data.json", "w") as out_file:
        json.dump(data, out_file, indent=4)


def read_data(name):
    with open("/Users/williamquinn/Desktop/PMT_Project/" + name + "_data.json", "r") as out_file:
        data = json.load(out_file)
    return data


def plot_charge(gamma_data, electron_data):
    plt.figure(figsize=figsize)
    freq, bin_edges = np.histogram(gamma_data["charge"], range=(0, 60), bins=240)
    width = bin_edges[2] - bin_edges[1]
    bin_centres = bin_edges[:-1] + width/2
    plt.bar(bin_centres, freq, width=width, alpha=0.5, label='gamma')
    freq, bin_edges = np.histogram(electron_data["charge"], range=(0, 60), bins=240)
    width = bin_edges[2] - bin_edges[1]
    bin_centres = bin_edges[:-1] + width / 2
    plt.bar(bin_centres, freq, width=width, alpha=0.5, label='electron')

    plt.xlabel("Charge /pC")
    plt.legend(loc='best')
    # plt.yscale('log')
    plt.xlim(0, 60)
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/PMT_Project/gamma_charge_comp.pdf")


def main():
    # gamma_file = "/Users/williamquinn/Desktop/PMT_Project/A1000_t1211_gamma.xml"
    # electron_file = "/Users/williamquinn/Desktop/PMT_Project/A1000_t1524.xml"
    # parse_data(gamma_file, "gamma")
    # parse_data(electron_file, "electron")
    gamma_data = read_data("gamma")
    result = root_fit_gamma(gamma_data, plot=True)
    electron_data = read_data("electron")
    root_fit_ce(electron_data, result, plot=True)
    # root_fit_both(electron_data, result, plot=True)
    # plot_charge(gamma_data=gamma_data, electron_data=electron_data)
    # print(">>> Got data")
    exit(1)


if __name__ == "__main__":
    main()
