import sys

import numpy as np
import tabulate

sys.path.insert(1, '../..')
from pmt_he_study.models import *


def print_cov_matrix(M, name):
    new_new_M = []
    mc = ['P0', 'P1', 'L']
    headers = [name, 'P0', 'P1', 'L']
    new_M = np.insert(M, 0, np.array(['1', '0', '2']), axis=1)
    for i in range(len(new_M)):
        new_new_M.append([])
        for j in range(len(new_M[i])):
            if j == 0:
                new_new_M[i].append(mc[i])
            else:
                new_new_M[i].append(str(new_M[i][j]))
    print(tabulate.tabulate(new_new_M, headers=headers, tablefmt="fancy_grid", numalign="right"))
    # print(tabulate.tabulate(new_new_M, headers=[name, 'm', 'c'], tablefmt="latex", numalign="right"))


def plot(model, dates, ap_charges, ap_charges_err, av_charges, av_charges_err, name):
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

    av_charge_0 = np.array(av_charges[0])
    av_charge_1 = np.array(av_charges[1])
    av_charge_err_0 = np.array(av_charges_err[0])
    av_charge_err_1 = np.array(av_charges_err[1])

    l = 5.91
    l_err = l/100
    y0 = (1 / 0.054) * ap_charge_0 / av_charge_0 / l
    y0_err = y0 * np.sqrt((ap_charge_err_0 / ap_charge_0) ** 2 + (av_charge_err_0 / av_charge_0) ** 2 + (l_err/l)**2)

    y1 = (1 / 0.054) * ap_charge_1 / av_charge_1 / l
    y1_err = y1 * np.sqrt((ap_charge_err_1 / ap_charge_1) ** 2 + (av_charge_err_1 / av_charge_1) ** 2 + (l_err/l)**2)

    x = np.array(date_0[mid + 1:]) - date_0[mid + 1:][0]
    y = y0[mid + 1:]
    yerr = y0_err[mid + 1:]
    pars, errs, chi, ndof = root_fit(x, y, yerr, model)
    # d_inv = np.diag(1 / np.sqrt(np.diag(pcov)))
    # pcor = d_inv @ pcov @ d_inv
    # print_cov_matrix(pcor, model.name + name)

    offset = pars[1]

    fig1 = plt.figure(num=None, figsize=(5, 4), dpi=80, facecolor='w', edgecolor='k')
    frame1 = fig1.add_axes((.15, .32, .8, .6))
    plt.errorbar(date_0[:start + 1], y0[:start + 1] -offset, zorder=0, yerr=y0_err[:start + 1], fmt="C0s",
                 label="Atmospheric He", markersize=1, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
    plt.errorbar(date_0[start + 1:mid + 1], y0[start + 1:mid + 1] -offset, zorder=0, yerr=y0_err[start + 1:mid + 1],
                 fmt="C1s", label="1% He", markersize=1, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
    plt.errorbar(date_0[mid + 1:], y0[mid + 1:] -offset, zorder=0, yerr=y0_err[mid + 1:], fmt="C2s", label="10% He",
                 markersize=1, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
    plt.errorbar(date_1, y1 -offset, zorder=0, yerr=y1_err, fmt="C3o", label="Control",
                 markersize=1, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
    plt.plot(date_0[mid + 1:], model.func(x, pars) -offset, 'k-', label='Model', zorder=10, linewidth=line_width)

    handles, labels = plt.gca().get_legend_handles_labels()
    n_sf = get_n_sfs(pars, errs)
    strings = [r'$P_0 =$ {:.' + str(n_sf[0]) + 'e} ± {:.0e}', r'$P_1 =$ {:.' + str(n_sf[1]) + 'f} ± {:.' + str(n_sf[1]) + 'f}',
               '$L =$ {:.0f} ± {:.0f} days']
    patch = patches.Patch(color='white', label=strings[0].format(pars[0], errs[0]))
    patch_1 = patches.Patch(color='white', label=strings[1].format(pars[1], errs[1]))
    patch_2 = patches.Patch(color='white',
                            label=strings[2].format(pars[2] / (3600 * 24), errs[2] / (3600 * 24)))
    patch_3 = patches.Patch(color='white', label=r'$\chi^2/N_{DoF}}$ = ' + '{:.2f}/{}'.format(chi*ndof, ndof))
    handles.extend([patch, patch_1, patch_2, patch_3])

    plt.ylabel(r'$p_i$ /Pa')
    plt.xlim(-30, 420)
    plt.legend(handles=handles, loc='upper left')

    frame2 = fig1.add_axes((.15, .1, .8, .2))
    plt.xlabel("Days from 1% Helium Onset")
    plt.axhline(0, ls='--', color='black')
    plt.ylabel("(model-data)/model")
    plt.xlim(-30, 420)
    plt.errorbar(date_0[mid + 1:],
                 (model.func(x, pars) - y0[mid + 1:]) / model.func(x, pars),
                 yerr=y0_err[mid + 1:] / model.func(date_0[mid + 1:], pars), fmt="k.",
                 markersize=1, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/PMT_Project/" + model.name + "_pi_vs_time_charge_" + name + ".pdf")
    plt.close()


def plot_ap_charge(model, dates, ap_charges, ap_charges_err, name):
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
    pars, errs, chi, ndof = root_fit(x, y, yerr, model)

    fig1 = plt.figure(num=None, figsize=(5, 5), dpi=80, facecolor='w', edgecolor='k')
    frame1 = fig1.add_axes((.15, .32, .8, .6))
    frame1.set_xticklabels([])
    plt.errorbar(date_0[:start + 1], ap_charge_0[:start + 1], zorder=0,
                 yerr=ap_charge_err_0[:start + 1], fmt="C0s", label="Atmospheric He",
                 markersize=1, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
    plt.errorbar(date_0[start + 1:mid + 1], ap_charge_0[start + 1:mid + 1], zorder=0,
                 yerr=ap_charge_err_0[start + 1:mid + 1], fmt="C1s", label="1% He",
                 markersize=1, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
    plt.errorbar(date_0[mid + 1:], ap_charge_0[mid + 1:], zorder=0,
                 yerr=ap_charge_err_0[mid + 1:], fmt="C2s", label="10% He",
                 markersize=1, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
    plt.errorbar(date_1, ap_charge_1, zorder=0,
                 yerr=ap_charge_err_1, fmt="C3o", label="Control",
                 markersize=1, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
    plt.plot(date_0[mid + 1:], model.func(x, pars), 'k-', label='Model', zorder=10, linewidth=line_width)
    # plt.axvline(date_0[start], 0, 100, ls='--', color='k')
    # plt.axvline(date_0[mid], 0, 100, ls='--', color='k')

    handles, labels = plt.gca().get_legend_handles_labels()

    n_sf = [0, 0, 0]
    for i in range(3):
        while errs[i] / pars[i] < 1:
            if errs[i] / pars[i] * pow(10, n_sf[i]) >= 1:
                break
            else:
                n_sf[i] += 1
    strings = [r'$P_0 =$ {:.' + str(n_sf[0]) + 'e} ± {:.0e}', r'$P_1 =$ {:.' + str(n_sf[1]) + 'e} ± {:.0e}',
               '$L =$ {:.0f} ± {:.0f} days']
    patch = patches.Patch(color='white', label=strings[0].format(pars[0], errs[0]))
    patch_1 = patches.Patch(color='white', label=strings[1].format(pars[1], errs[1]))
    patch_2 = patches.Patch(color='white',
                            label=strings[2].format(pars[2] / (3600 * 24), errs[2] / (3600 * 24)))
    patch_3 = patches.Patch(color='white', label=r'$\chi^2_R =$ {:.2f}'.format(chi))
    handles.extend([patch, patch_1, patch_2, patch_3])

    plt.ylabel("Average Charge /pC")
    plt.title("Average After-pulse Region Charge")
    plt.xlim(-30, 420)
    plt.legend(handles=handles, loc='upper left')

    frame2 = fig1.add_axes((.15, .1, .8, .2))
    plt.xlabel("Days from 1% Helium Onset")
    plt.axhline(0, ls='--', color='black')
    plt.ylabel("(model-data)/model")
    plt.xlim(-30, 420)
    plt.errorbar(date_0[mid + 1:],
                 (model.func(x, pars) - ap_charge_0[mid + 1:]) / model.func(x, pars),
                 yerr=ap_charge_err_0[mid + 1:]/model.func(date_0[mid + 1:], pars), fmt="k.",
                 markersize=1, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/PMT_Project/" + model.name + "ap_charge_vs_time_" + name + ".pdf")
    plt.close()


def plot_charge_ratio(model, dates, ratios, ratios_err, name: str):
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
    pars, errs, chi, ndof = root_fit(x, y, yerr, model)

    fig1 = plt.figure(num=None, figsize=(5,5), dpi=80, facecolor='w', edgecolor='k')
    frame1 = fig1.add_axes((.15, .32, .8, .6))
    frame1.set_xticklabels([])
    plt.errorbar(date_0[:start + 1], ratio_0[:start + 1], zorder=0,
                 yerr=ratio_err_0[:start + 1], fmt="C0s", label="Atmospheric He",
                 markersize=1, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
    plt.errorbar(date_0[start + 1:mid + 1], ratio_0[start + 1:mid + 1], zorder=0,
                 yerr=ratio_err_0[start + 1:mid + 1], fmt="C1s", label="1% He",
                 markersize=1, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
    plt.errorbar(date_0[mid + 1:], ratio_0[mid + 1:], zorder=0,
                 yerr=ratio_err_0[mid + 1:], fmt="C2s", label="10% He",
                 markersize=1, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
    plt.errorbar(date_1, ratio_1, zorder=0,
                 yerr=ratio_err_1, fmt="C3o", label="Control",
                 markersize=1, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
    plt.plot(date_0[mid + 1:], model.func(x, pars), 'k-', label='Model', zorder=10, linewidth=line_width)
    # plt.axvline(date_0[start], 0, 100, ls='--', color='k')
    # plt.axvline(date_0[mid], 0, 100, ls='--', color='k')

    handles, labels = plt.gca().get_legend_handles_labels()

    n_sf = [0, 0, 0]
    for i in range(3):
        while errs[i] / pars[i] < 1:
            if errs[i] / pars[i] * pow(10, n_sf[i]) >= 1:
                break
            else:
                n_sf[i] += 1
    strings = [r'$P_2 =$ {:.' + str(n_sf[0]) + 'e} ± {:.0e}', r'$P_3 =$ {:.' + str(n_sf[1]) + 'e} ± {:.0e}',
               '$L =$ {:.0f} ± {:.0f} days']
    patch = patches.Patch(color='white', label=strings[0].format(pars[0], errs[0]))
    patch_1 = patches.Patch(color='white', label=strings[1].format(pars[1], errs[1]))
    patch_2 = patches.Patch(color='white',
                            label=strings[2].format(pars[2] / (3600 * 24), errs[2] / (3600 * 24)))
    patch_3 = patches.Patch(color='white', label=r'$\chi^2_R =$ {:.2f}'.format(chi))
    handles.extend([patch, patch_1, patch_2, patch_3])

    plt.ylabel("Ratio Average")
    plt.title("After-pulse Region Charge vs Pulse Charge Ratio")
    plt.xlim(-30, 420)
    plt.legend(handles=handles, loc='upper left')

    frame2 = fig1.add_axes((.15, .1, .8, .2))
    plt.xlabel("Days from 10% He Onset")
    plt.axhline(0, ls='--', color='black')
    plt.ylabel("(model-data)/model")
    plt.xlim(-30, 420)
    plt.errorbar(date_0[mid + 1:],
                 (model.func(x, pars) - ratio_0[mid + 1:]) / model.func(x, pars),
                 yerr=ratio_err_0[mid + 1:] / model.func(date_0[mid + 1:], pars), fmt="k.",
                 markersize=1, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
    plt.savefig("/Users/williamquinn/Desktop/PMT_Project/" + model.name + "_ratio_vs_time_" + name + ".pdf")


    # extrapolate(model, pars, 1, name)

    '''plt.figure(figsize=figsize)
    plt.plot(x, y)
    plt.xlim(0, x[-1])
    plt.title("Charge Ratio Model Projection")
    plt.xlabel("Days from 10% He Onset")
    plt.ylabel("Ratio")
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/PMT_Project/" + model.name + "_ratio_projection_" + name + ".pdf")
    plt.close()
'''


def fit_ratios(hist):
    freq = []
    bin_centres = []
    for i_bin in range(hist.GetNbinsX()):
        if hist.GetBinContent(i_bin + 1) > 0:
            freq.append(hist.GetBinContent(i_bin + 1))
            bin_centres.append(hist.GetBinCenter(i_bin + 1))
    width = hist.GetBinWidth(1)
    freq = np.array(freq)
    bin_centres = np.array(bin_centres)
    # plt.bar(bin_centres, freq, width=width, alpha=0.5)
    # plt.axvline(hist.GetMean(), ls='--')
    # plt.yscale('log')
    # plt.show()

    mean = hist.GetMean()
    mean_err = hist.GetMeanError()

    '''try:
        popt, pcov = curve_fit(f=lorentzian, xdata=bin_centres, ydata=freq, p0=[100, 0, 0.1],
                               bounds=[[0, -1, 0], [1e6, 1, 10]], maxfev=1000000)
        print("Fitting Lorenztian to ratio...")
    except RuntimeError:
        popt, pcov = curve_fit(f=gaussian_noh, xdata=bin_centres, ydata=freq, p0=[100, 0, 0.1],
                               bounds=[[0, -1, 0], [1e6, 1, 10]], maxfev=10000)
        print("Fitting Gaussian to ratio...")
    perr = np.sqrt(np.diag(pcov))
    par_names = ['A', 'loc', 'scale']
    chi_2, ndof = chi2(freq, np.sqrt(freq), lorentzian(bin_centres, *popt), len(popt))
    for i in range(len(popt)):
        print("{}: {} ± {}".format(par_names[i], popt[i], perr[i]))
    print("chi2/Ndof:", "{:.2f}/{}".format(chi_2, ndof))'''

    return mean, mean_err


def main():
    topology = [2, 1]
    pmt_array = PMT_Array(topology, "summary")
    pmt_array.set_pmt_id("GAO607", 0)
    pmt_array.set_pmt_id("GAO612", 1)

    filenames_txt = "/Users/williamquinn/Desktop/PMT_Project/data/1400V/filenames.txt"
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

    av_charge = [[], []]
    av_charge_err = [[], []]

    ap_charge_err = [[], []]
    ap_ratio_err = [[], []]
    he_ap_charge_err = [[], []]
    he_ap_ratio_err = [[], []]

    for i_file in range(filenames.size):
        filename = filenames[i_file][0].decode("utf-8")
        date = filename.split("_")[0]
        voltage = int(filename.split("_")[1].split("A")[1])

        file = ROOT.TFile("/Users/williamquinn/Desktop/PMT_Project/data/1400V/" + filename, "READ")
        file.cd()

        for i_om in range(2):
            charge_hist = file.Get(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() +
                                   "_charge_spectrum_" + str(voltage) + "V")
            ap_charge_hist = file.Get(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() +
                                      "_ap_charge_spectrum_" + str(voltage) + "V")
            he_ap_charge_hist = file.Get(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() +
                                         "_he_ap_charge_spectrum_" + str(voltage) + "V")
            ap_charge_charge_hist = file.Get(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() +
                                             "_ap_charge_charge_spectrum_" + str(voltage) + "V")
            he_ap_charge_charge_hist = file.Get(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() +
                                                "_he_ap_charge_charge_spectrum_" + str(voltage) + "V")

            try:
                charge_hist.GetEntries()
                ap_charge_hist.GetEntries()
                he_ap_charge_hist.GetEntries()
                ap_charge_charge_hist.GetEntries()
                he_ap_charge_charge_hist.GetEntries()
            except:
                continue

            if charge_hist.GetMean() < 20:
                continue

            dates[i_om].append(int(date))

            ap_charge[i_om].append(ap_charge_hist.GetMean())
            ap_charge_err[i_om].append(ap_charge_hist.GetMeanError())
            he_ap_charge[i_om].append(he_ap_charge_hist.GetMean())
            he_ap_charge_err[i_om].append(he_ap_charge_hist.GetMeanError())
            av_charge[i_om].append(charge_hist.GetMean())
            av_charge_err[i_om].append(charge_hist.GetMeanError())

            ratio_mean, ratio_mean_err = fit_ratios(ap_charge_charge_hist)
            he_ratio_mean, he_ratio_mean_err = fit_ratios(he_ap_charge_charge_hist)
            ap_ratio[i_om].append(ratio_mean)
            he_ap_ratio[i_om].append(he_ratio_mean)
            ap_ratio_err[i_om].append(ratio_mean_err)
            he_ap_ratio_err[i_om].append(he_ratio_mean_err)

            del ap_charge_hist
            del he_ap_charge_hist
            del ap_charge_charge_hist
            del he_ap_charge_charge_hist

    model = Model()
    # plot_ap_charge(model, dates, ap_charge, ap_charge_err, "")
    # plot_ap_charge(model, dates, he_ap_charge, he_ap_charge_err, "he")
    # plot_charge_ratio(model, dates, ap_ratio, ap_ratio_err, "")
    # plot_charge_ratio(model, dates, he_ap_ratio, he_ap_ratio_err, "he")
    # plot(model, dates, ap_charge, ap_charge_err, av_charge, av_charge_err, "")
    plot(model, dates, he_ap_charge, he_ap_charge_err, av_charge, av_charge_err, "he")

    model = Model_0()
    # plot_ap_charge(model, dates, ap_charge, ap_charge_err, "")
    # plot_ap_charge(model, dates, he_ap_charge, he_ap_charge_err, "he")
    # plot_charge_ratio(model, dates, ap_ratio, ap_ratio_err, "")
    # plot_charge_ratio(model, dates, he_ap_ratio, he_ap_ratio_err, "he")
    # plot(model, dates, ap_charge, ap_charge_err, av_charge, av_charge_err, "")
    plot(model, dates, he_ap_charge, he_ap_charge_err, av_charge, av_charge_err, "he")


if __name__ == "__main__":
    main()
