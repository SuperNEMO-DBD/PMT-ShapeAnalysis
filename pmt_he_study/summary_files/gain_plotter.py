import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(1, '../..')
from pmt_he_study.models import *

e = 1.603E-19
names = ['Exposed', 'Control']


def sort_dict(the_dict: dict):
    new_w0 = [i for i in range(-5, 59)]
    new_dict = the_dict.copy()

    for volt in [1000, 1400]:
        if new_dict[volt] is not None:
            for ch in [0, 1]:
                for i in range(len(new_w0)):
                    if new_dict[volt][ch]["w"][i] != new_w0[i]:
                        new_dict[volt][ch]["w"] = np.insert(new_dict[volt][ch]["w"], i, new_w0[i])
                        new_dict[volt][ch]["v"] = np.insert(new_dict[volt][ch]["v"], i, -9999)
                        new_dict[volt][ch]["e"] = np.insert(new_dict[volt][ch]["e"], i, -9999)
                    else:
                        if new_dict[volt][ch]["e"][i] < new_dict[volt][ch]["v"][i] * 0.001:
                            new_dict[volt][ch]["e"][i] = new_dict[volt][ch]["v"][i] * 0.001
    return new_dict


def write_to_file(filename: str, av_charge, mev_charge, gain, res):
    w = [i for i in range(-5, 59)]
    av_charge_d = {
        "name": "av_charge",
        1000: {
            0: {
                "w": av_charge[0][0][0],
                "v": av_charge[0][1][0],
                "e": av_charge[0][2][0]
            },
            1: {
                "w": av_charge[0][0][1],
                "v": av_charge[0][1][1],
                "e": av_charge[0][2][1]
            }
        },
        1400: {
            0: {
                "w": av_charge[1][0][0],
                "v": av_charge[1][1][0],
                "e": av_charge[1][2][0]
            },
            1: {
                "w": av_charge[1][0][1],
                "v": av_charge[1][1][1],
                "e": av_charge[1][2][1]
            }
        }
    }
    mev_charge_d = {
        "name": "mev_charge",
        1000: {
            0: {
                "w": mev_charge[0][0][0],
                "v": mev_charge[0][1][0],
                "e": mev_charge[0][2][0]
            },
            1: {
                "w": mev_charge[0][0][1],
                "v": mev_charge[0][1][1],
                "e": mev_charge[0][2][1]
            }
        },
        1400: {
            0: {
                "w": mev_charge[1][0][0],
                "v": mev_charge[1][1][0],
                "e": mev_charge[1][2][0]
            },
            1: {
                "w": mev_charge[1][0][1],
                "v": mev_charge[1][1][1],
                "e": mev_charge[1][2][1]
            }
        }
    }
    gain_d = {
        "name": "gain",
        1000: {
            0: {
                "w": gain[0][0],
                "v": gain[1][0],
                "e": gain[2][0]
            },
            1: {
                "w": gain[0][1],
                "v": gain[1][1],
                "e": gain[2][1]
            }
        },
        1400: None
    }
    res_d = {
        "name": "res",
        1000: {
            0: {
                "w": res[0][0],
                "v": res[1][0],
                "e": res[2][0]
            },
            1: {
                "w": res[0][1],
                "v": res[1][1],
                "e": res[2][1]
            }
        },
        1400: None
    }

    new_av_charge = sort_dict(av_charge_d)
    new_mev_charge = sort_dict(mev_charge_d)
    new_gain = sort_dict(gain_d)
    new_res = sort_dict(res_d)

    file = open(filename, "w")
    column_headers = ['week',
                      'Ch0 MeV Charge at 1000V', 'err',
                      'Ch1 MeV Charge at 1000V', 'err',
                      'Ch0 MeV Charge at 1400V', 'err',
                      'Ch1 MeV Charge at 1400V', 'err',
                      'Ch0 Average Charge at 1000V', 'err',
                      'Ch1 Average Charge at 1000V', 'err',
                      'Ch0 Average Charge at 1400V', 'err',
                      'Ch1 Average Charge at 1400V', 'err',
                      'Ch0 gain at 1000V', 'err',
                      'Ch1 gain at 1000V', 'err',
                      'Ch0 gain at 1400V', 'err',
                      'Ch1 gain at 1400V', 'err',
                      'Ch0 Resolution at 1000V', 'err',
                      'Ch1 Resolution at 1000V', 'err'
                      ]
    header_string = ''
    for col in column_headers:
        header_string += col + ","
    header_string = header_string[:-1] + '\n'
    file.write(header_string)

    for i in range(len(w)):
        g_14_0 = new_gain[1000][0]["v"][i] * new_mev_charge[1400][0]["v"][i]/new_mev_charge[1000][0]["v"][i]
        g_14_1 = new_gain[1000][1]["v"][i] * new_mev_charge[1400][1]["v"][i]/new_mev_charge[1000][1]["v"][i]
        file.write("{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n".format(w[i],
                                                                                                             new_mev_charge[1000][0]["v"][i], new_mev_charge[1000][0]["e"][i],
                                                                                                             new_mev_charge[1000][1]["v"][i], new_mev_charge[1000][1]["e"][i],
                                                                                                             new_mev_charge[1400][0]["v"][i], new_mev_charge[1400][0]["e"][i],
                                                                                                             new_mev_charge[1400][1]["v"][i], new_mev_charge[1400][1]["e"][i],
                                                                                                             new_av_charge[1000][0]["v"][i], new_av_charge[1000][0]["e"][i],
                                                                                                             new_av_charge[1000][1]["v"][i], new_av_charge[1000][1]["e"][i],
                                                                                                             new_av_charge[1400][0]["v"][i], new_av_charge[1400][0]["e"][i],
                                                                                                             new_av_charge[1400][1]["v"][i], new_av_charge[1400][1]["e"][i],
                                                                                                             new_gain[1000][0]["v"][i], new_gain[1000][0]["e"][i],
                                                                                                             new_gain[1000][1]["v"][i], new_gain[1000][1]["e"][i],
                                                                                                             new_gain[1000][0]["v"][i] * new_mev_charge[1400][0]["v"][i]/new_mev_charge[1000][0]["v"][i],
                                                                                                             new_gain[1000][0]["v"][i] * new_mev_charge[1400][0]["v"][i]/new_mev_charge[1000][0]["v"][i] * np.sqrt((new_gain[1000][0]["e"][i]/new_gain[1000][0]["v"][i])**2 + (new_mev_charge[1000][0]["e"][i]/new_mev_charge[1000][0]["v"][i])**2 + (new_mev_charge[1400][0]["e"][i]/new_mev_charge[1400][0]["v"][i])**2),
                                                                                                             new_gain[1000][1]["v"][i] * new_mev_charge[1400][1]["v"][i]/new_mev_charge[1000][1]["v"][i],
                                                                                                             new_gain[1000][0]["v"][i] * new_mev_charge[1400][0]["v"][i]/new_mev_charge[1000][0]["v"][i] * np.sqrt((new_gain[1000][1]["e"][i]/new_gain[1000][1]["v"][i])**2 + (new_mev_charge[1000][1]["e"][i]/new_mev_charge[1000][1]["v"][i])**2 + (new_mev_charge[1400][1]["e"][i]/new_mev_charge[1400][1]["v"][i])**2),
                                                                                                             new_res[1000][0]["v"][i], new_res[1000][0]["e"][i],
                                                                                                             new_res[1000][1]["v"][i], new_res[1000][1]["e"][i]
                                                                                                             ))
    file.close()


def plot_gain(data, voltage):
    charge_cut = 0
    if voltage == 1400:
        charge_cut = 200
    else:
        charge_cut = 20

    x, g, g_err = [[], []], [[], []], [[], []]
    for k in range(2):
        x[k] = [val for index, val in enumerate(data["dates"][k]) if
                data["chi2s"][k][index] < 10 and data["mev_charge"][k][index] > charge_cut]
        g[k] = [val for index, val in enumerate(data["gain"][k]) if
                data["chi2s"][k][index] < 10 and data["mev_charge"][k][index] > charge_cut]
        g_err[k] = [val for index, val in enumerate(data["gain_err"][k]) if
                    data["chi2s"][k][index] < 10 and data["mev_charge"][k][index] > charge_cut]

    plt.figure(figsize=figsize, facecolor='white')
    for i in range(2):
        plt.errorbar(x[i], g[i], yerr=g_err[i], fmt='.', label=names[i], markersize=marker_size, capsize=cap_size,
                     linewidth=line_width, capthick=cap_thick)

    plt.legend(loc='best')
    plt.xlabel("Days from onset of 1% Helium")
    plt.ylabel("Gain ")
    plt.title("Gain at {}V".format(voltage))
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/PMT_Project/gain_{}V.pdf".format(voltage))

    w = [[i for i in range(-5, 59)], [i for i in range(-5, 59)]]
    av_g = [[np.array([]) for i in range(-5, 59)], [np.array([]) for i in range(-5, 59)]]
    av_g_err = [[np.array([]) for i in range(-5, 59)], [np.array([]) for i in range(-5, 59)]]

    for i in range(2):
        for index, val in enumerate(x[i]):
            week = val // 7
            pos = week + 5
            av_g[i][pos] = np.append(av_g[i][pos], g[i][index])
            av_g_err[i][pos] = np.append(av_g_err[i][pos], g_err[i][index])

    w_ = [[], []]
    av_gain = [[], []]
    av_gain_err = [[], []]
    for i in range(2):
        w_[i] = np.array([val for index, val in enumerate(w[i]) if len(av_g[i][index]) > 0])
        av_gain[i] = np.array([np.average(val, weights=av_g_err[i][index]) for index, val in enumerate(av_g[i]) if len(val) > 0])
        av_gain_err[i] = np.array([np.sqrt(np.average((val - np.average(val, weights=av_g_err[i][index]))**2, weights=av_g_err[i][index])/len(val)) for index, val in enumerate(av_g[i]) if len(val) > 0])

    plt.figure(figsize=figsize, facecolor='white')
    for i in range(2):
        plt.errorbar(x[i], g[i], yerr=g_err[i], fmt=".", label=names[i], markersize=marker_size, capsize=cap_size,
                     linewidth=line_width, capthick=cap_thick)
        plt.errorbar(np.array(w_[i]) * 7 + 3, av_gain[i], yerr=av_gain_err[i], fmt=".", label=names[i] + " week",
                     markersize=marker_size, capsize=cap_size, linewidth=line_width, capthick=cap_thick)

    plt.legend(loc='best')
    plt.xlabel("Days from onset of 1% Helium")
    plt.ylabel("Gain ")
    plt.title("Gain at {}V".format(voltage))
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/PMT_Project/gain_{}V_week.pdf".format(voltage))

    return w_, av_gain, av_gain_err


def plot_res(data):
    charge_cut = 20
    x, res, res_err = [[], []], [[], []], [[], []]
    for k in range(2):
        x[k] = [val for index, val in enumerate(data["dates"][k]) if
                data["chi2s"][k][index] < 10 and data["mev_charge"][k][index] > charge_cut]
        res[k] = [val for index, val in enumerate(data["res"][k]) if
                  data["chi2s"][k][index] < 10 and data["mev_charge"][k][index] > charge_cut]
        res_err[k] = [val for index, val in enumerate(data["res_err"][k]) if
                      data["chi2s"][k][index] < 10 and data["mev_charge"][k][index] > charge_cut]

    plt.figure(figsize=figsize, facecolor='white')
    for i in range(2):
        plt.errorbar(x[i], res[i], yerr=res_err[i], fmt='.', label=names[i], markersize=marker_size,
                     capsize=cap_size, linewidth=line_width, capthick=cap_thick)

    plt.legend(loc='best')
    plt.xlabel("Days from onset of 1% Helium")
    plt.ylabel("Resoultion")
    plt.title("Resolution at 1MeV")
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/PMT_Project/res.pdf")

    # av = [[0 for i in range(-5, 59)], [0 for i in range(-5, 59)]]
    w = [[i for i in range(-5, 59)], [i for i in range(-5, 59)]]
    av_r = [[np.array([]) for i in range(-5, 59)], [np.array([]) for i in range(-5, 59)]]
    av_r_err = [[np.array([]) for i in range(-5, 59)], [np.array([]) for i in range(-5, 59)]]

    for i in range(2):
        for index, val in enumerate(x[i]):
            week = val // 7
            pos = week + 5
            av_r[i][pos] = np.append(av_r[i][pos], res[i][index])
            av_r_err[i][pos] = np.append(av_r_err[i][pos], res_err[i][index])

    w_ = [[], []]
    av_res = [[], []]
    av_res_err = [[], []]
    for i in range(2):
        w_[i] = np.array([val for index, val in enumerate(w[i]) if len(av_r[i][index]) > 0])
        av_res[i] = np.array([np.average(val, weights=av_r_err[i][index]) for index, val in enumerate(av_r[i]) if len(val) > 0])
        av_res_err[i] = np.array([np.sqrt(np.average((val - np.average(val, weights=av_r_err[i][index]))**2, weights=av_r_err[i][index]) / len(val)) for index, val in enumerate(av_r[i]) if len(val) > 0])

    plt.figure(figsize=figsize, facecolor='white')
    for i in range(2):
        plt.errorbar(x[i], res[i], yerr=res_err[i], fmt=".", label=names[i], markersize=marker_size, capsize=cap_size,
                     linewidth=line_width, capthick=cap_thick)
        plt.errorbar(w_[i] * 7 + 3, av_res[i], yerr=av_res_err[i], fmt=".", label=names[i] + " week",
                     markersize=marker_size, capsize=cap_size, linewidth=line_width, capthick=cap_thick)

    plt.legend(loc='best')
    plt.xlabel("Days from onset of 1% Helium")
    plt.ylabel("Resoultion")
    plt.title("Resolution at 1MeV")
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/PMT_Project/res_week.pdf")

    return w_, av_res, av_res_err


def plot_mev_charge(data, voltage):

    charge_cut = 0
    if voltage == 1400:
        charge_cut = 200
    else:
        charge_cut = 20
    x, y, y_err = [[], []], [[], []], [[], []]
    for k in range(2):
        x[k] = [val for index, val in enumerate(data["dates"][k]) if data["chi2s"][k][index] < 10 and data["mev_charge"][k][index] > charge_cut]
        y[k] = [val for index, val in enumerate(data["mev_charge"][k]) if data["chi2s"][k][index] < 10 and val > charge_cut]
        y_err[k] = [val for index, val in enumerate(data["mev_charge_err"][k]) if data["chi2s"][k][index] < 10 and data["mev_charge"][k][index] > charge_cut]

    plt.figure(figsize=figsize, facecolor='white')
    for i in range(2):
        plt.errorbar(x[i], y[i], yerr=y_err[i], fmt='.', label=names[i], markersize=marker_size, capsize=cap_size,
                     linewidth=line_width, capthick=cap_thick)

    plt.legend(loc='best')
    plt.xlabel("Days from onset of 1% Helium")
    plt.ylabel("Charge at 1MeV /pC")
    plt.title("Charge at 1MeV vs Time at {}V".format(voltage))
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/PMT_Project/MeV_charge_vs_time_{}V.pdf".format(voltage))

    w = [[i for i in range(-5, 59)], [i for i in range(-5, 59)]]
    av_y = [[np.array([]) for i in range(-5, 59)], [np.array([]) for i in range(-5, 59)]]
    av_y_err = [[np.array([]) for i in range(-5, 59)], [np.array([]) for i in range(-5, 59)]]

    for i in range(2):
        for index, val in enumerate(x[i]):
            week = val // 7
            pos = week + 5
            av_y[i][pos] = np.append(av_y[i][pos], y[i][index])
            av_y_err[i][pos] = np.append(av_y_err[i][pos], y_err[i][index])

    w_ = [[], []]
    av_mev = [[], []]
    av_mev_err = [[], []]
    for i in range(2):
        w_[i] = np.array([val for index, val in enumerate(w[i]) if len(av_y[i][index]) > 0])
        av_mev[i] = np.array([np.average(val, weights=av_y_err[i][index]) for index, val in enumerate(av_y[i]) if len(val) > 0])
        av_mev_err[i] = np.array([np.sqrt(np.average((val - np.average(val, weights=av_y_err[i][index]))**2, weights=av_y_err[i][index])/len(val)) for index, val in enumerate(av_y[i]) if len(val) > 0])

    plt.figure(figsize=figsize, facecolor='white')
    for i in range(2):
        plt.errorbar(x[i], y[i], yerr=y_err[i], fmt=".", label=names[i], markersize=marker_size, capsize=cap_size,
                     linewidth=line_width, capthick=cap_thick)
        plt.errorbar(w_[i] * 7 + 3, av_mev[i], yerr=av_mev_err[i], fmt=".", label=names[i] + " week",
                     markersize=marker_size, capsize=cap_size, linewidth=line_width, capthick=cap_thick)

    plt.legend(loc='best')
    plt.xlabel("Days from onset of 1% Helium")
    plt.ylabel("Charge at 1MeV /pC")
    plt.title("Charge at 1MeV vs Time at {}V".format(voltage))
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/PMT_Project/MeV_charge_vs_time_{}V_week.pdf".format(voltage))

    return w_, av_mev, av_mev_err


def plot_av_charge(data, voltage):
    charge_cut = 0
    if voltage == 1400:
        charge_cut = 100
    else:
        charge_cut = 10

    x, y, y_err = [[], []], [[], []], [[], []]
    for k in range(2):
        x[k] = [val for index, val in enumerate(data["dates"][k]) if data["av_charge"][k][index] > charge_cut]
        y[k] = [val for index, val in enumerate(data["av_charge"][k]) if val > charge_cut]
        y_err[k] = [val for index, val in enumerate(data["av_charge_err"][k]) if data["av_charge"][k][index] > charge_cut]

    plt.figure(figsize=figsize)
    for i in range(2):
        plt.errorbar(x[i], y[i], yerr=y_err[i], fmt=".", label=names[i], markersize=marker_size, capsize=cap_size,
                     linewidth=line_width, capthick=cap_thick)
    plt.legend(loc='best')
    plt.xlabel("Days from onset of 1% Helium")
    plt.ylabel("Charge /pC")
    '''if voltage == 1400:
        plt.ylim(110, 150)'''
    plt.title("Average Pulse Charge at {}V".format(voltage))
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/PMT_Project/av_pulse_charge_vs_time_{}V.pdf".format(voltage))

    w = [[i for i in range(-5, 59)], [i for i in range(-5, 59)]]
    av_y = [[np.array([]) for i in range(-5, 59)], [np.array([]) for i in range(-5, 59)]]
    av_y_err = [[np.array([]) for i in range(-5, 59)], [np.array([]) for i in range(-5, 59)]]

    for i in range(2):
        for index, val in enumerate(x[i]):
            week = val // 7
            pos = week + 5
            av_y[i][pos] = np.append(av_y[i][pos], y[i][index])
            av_y_err[i][pos] = np.append(av_y_err[i][pos], y_err[i][index])

    w_ = [[], []]
    av_tot = [[], []]
    av_tot_err = [np.array([]), np.array([])]
    for i in range(2):
        w_[i] = np.array([val for index, val in enumerate(w[i]) if len(av_y[i][index]) > 0])
        av_tot[i] = np.array([np.average(val, weights=av_y_err[i][index]) for index, val in enumerate(av_y[i]) if len(val) > 0])
        av_tot_err[i] = np.array([np.sqrt(np.average((val - np.average(val, weights=av_y_err[i][index]))**2, weights=av_y_err[i][index])/len(val)) for index, val in enumerate(av_y[i]) if len(val) > 0])

    plt.figure(figsize=figsize)
    for i in range(2):
        plt.errorbar(x[i], y[i], yerr=y_err[i], fmt=".", label=names[i], markersize=marker_size, capsize=cap_size,
                     linewidth=line_width, capthick=cap_thick)
        plt.errorbar(w_[i] * 7 + 3, av_tot[i], yerr=av_tot_err[i], fmt=".", label=names[i] + " week",
                     markersize=marker_size, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
    plt.legend(loc='best')
    plt.xlabel("Days from onset of 1% Helium")
    plt.ylabel("Charge /pC")
    '''if voltage == 1400:
        plt.ylim(110, 150)'''
    plt.title("Average Pulse Charge at {}V".format(voltage))
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/PMT_Project/av_pulse_charge_vs_time_{}V_week.pdf".format(voltage))

    return w_, av_tot, av_tot_err


def read_1400V_files(pmt_array: PMT_Array, filenames: str):
    try:
        print(">>> Reading data from file: {}".format(filenames))
        date_file = open(filenames, 'r')
    except FileNotFoundError as fnf_error:
        print(fnf_error)
        raise Exception("Error opening data file {}".format(filenames))

    filenames = np.loadtxt(filenames, delimiter=',', dtype={
        'names': ['filename'],
        'formats': ['S100']}, unpack=True)

    dates = [[], []]
    av_charges = [[], []]
    av_charges_err = [[], []]
    mev_charges = [[], []]
    mev_charges_err = [[], []]
    chi2s = [[], []]
    for i_file in range(filenames.size):
        filename = filenames[i_file][0].decode("utf-8")
        date = filename.split("_")[0]
        voltage = int(filename.split("_")[1].split("A")[1])

        file = ROOT.TFile("/Users/williamquinn/Desktop/data/1400V/" + filename, "READ")
        file.cd()

        for i_om in range(2):
            charge_hist = file.Get(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() +
                                   "_charge_spectrum_" + str(voltage) + "V")
            try:
                charge_hist.GetEntries()
            except:
                continue
            mean = get_sat_bi_mean(charge_hist)
            charge_pars = do_bi_1400(charge_hist, mean)

            dates[i_om].append(int(date))

            av_charges[i_om].append(charge_hist.GetMean())
            av_charges_err[i_om].append(charge_hist.GetMeanError())

            mev_charges[i_om].append(charge_pars[0][0])
            mev_charges_err[i_om].append(charge_pars[0][1])

            chi2s[i_om].append(charge_pars[-1])

            del charge_hist

        file.Close()

    data = {
        "dates": [process_date(dates[0]), process_date(dates[1])],
        "av_charge": [np.array(av_charges[0]), np.array(av_charges[1])],
        "av_charge_err": [np.array(av_charges_err[0]), np.array(av_charges_err[1])],
        "mev_charge": [np.array(mev_charges[0]), np.array(mev_charges[1])],
        "mev_charge_err": [np.array(mev_charges_err[0]), np.array(mev_charges_err[1])],
        "chi2s": [np.array(chi2s[0]), np.array(chi2s[1])]
    }

    return data


def read_1000V_files(pmt_array: PMT_Array, filenames: str):

    try:
        print(">>> Reading data from file: {}".format(filenames))
        date_file = open(filenames, 'r')
    except FileNotFoundError as fnf_error:
        print(fnf_error)
        raise Exception("Error opening data file {}".format(filenames))

    filenames = np.loadtxt(filenames, delimiter=',', dtype={
        'names': ['filename'],
        'formats': ['S100']}, unpack=True)

    dates = [[], []]
    mev_charges = [[], []]
    mev_charges_err = [[], []]
    av_charges = [[], []]
    av_charges_err = [[], []]
    res = [[], []]
    res_err = [[], []]
    gains = [[], []]
    gains_err = [[], []]
    chi2s = [[], []]

    for i_file in range(filenames.size):
        filename = filenames[i_file][0].decode("utf-8")
        date = filename.split("_")[0]
        voltage = int(filename.split("_")[1].split("A")[1])

        file = ROOT.TFile("/Users/williamquinn/Desktop/data/1000V/" + filename, "READ")
        file.cd()

        for i_om in range(2):
            charge_hist = file.Get(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() +
                                   "_charge_spectrum_" + str(voltage) + "V")

            try:
                charge_hist.GetEntries()
            except:
                continue
            mean = get_sat_bi_mean(charge_hist)

            try:
                charge_pars = do_bi_1000(charge_hist, mean)
            except ValueError:
                print("")
                continue

            dates[i_om].append(int(date))

            mev_charges[i_om].append(charge_pars[0][0])
            mev_charges_err[i_om].append(charge_pars[0][1])

            av_charges[i_om].append(charge_hist.GetMean())
            av_charges_err[i_om].append(charge_hist.GetMeanError())

            R = charge_pars[1][0] / charge_pars[0][0]
            G = (R ** 2) * charge_pars[0][0] * 1e-12 / e
            gains[i_om].append(G)
            gains_err[i_om].append(G * np.sqrt( (charge_pars[0][1]/charge_pars[0][0])**2 + 2*(charge_pars[1][1]/charge_pars[1][0])**2))

            res[i_om].append(R)
            res_err[i_om].append(R*np.sqrt((charge_pars[0][1]/charge_pars[0][0])**2 + (charge_pars[1][1]/charge_pars[1][0])**2))

            chi2s[i_om].append(charge_pars[-1])

            del charge_hist

    data = {
        "dates": [process_date(dates[0]), process_date(dates[1])],
        "av_charge": [np.array(av_charges[0]), np.array(av_charges[1])],
        "av_charge_err": [np.array(av_charges_err[0]), np.array(av_charges_err[1])],
        "mev_charge": [np.array(mev_charges[0]), np.array(mev_charges[1])],
        "mev_charge_err": [np.array(mev_charges_err[0]), np.array(mev_charges_err[1])],
        "res": [np.array(res[0]), np.array(res[1])],
        "res_err": [np.array(res_err[0]), np.array(res_err[1])],
        "gain": [np.array(gains[0]), np.array(gains[1])],
        "gain_err": [np.array(gains_err[0]), np.array(gains_err[1])],
        "chi2s": [np.array(chi2s[0]), np.array(chi2s[1])]
    }

    return data


def main():
    topology = [2, 1]
    pmt_array = PMT_Array(topology, "summary")
    pmt_array.set_pmt_id("GAO612", 1)
    pmt_array.set_pmt_id("GAO607", 0)

    filenames_txt = "/Users/williamquinn/Desktop/data/1400V/filenames.txt"
    data_1400V = read_1400V_files(pmt_array, filenames_txt)

    av_charge_14kV = plot_av_charge(data_1400V, 1400)
    mev_charge_14kV = plot_mev_charge(data_1400V, 1400)

    filenames_txt = "/Users/williamquinn/Desktop/data/1000V/filenames.txt"
    data_1000V = read_1000V_files(pmt_array, filenames_txt)

    av_charge_1kV = plot_av_charge(data_1000V, 1000)
    mev_charge_1kV = plot_mev_charge(data_1000V, 1000)
    gain_1kV = plot_gain(data_1000V, 1000)
    res_1kV = plot_res(data_1000V)

    av_charge = [av_charge_1kV, av_charge_14kV]
    mev_charge = [mev_charge_1kV, mev_charge_14kV]

    output_file = "/Users/williamquinn/Desktop/PMT_Project/week_values.csv"
    write_to_file(output_file, av_charge, mev_charge, gain_1kV, res_1kV)


if __name__ == "__main__":
    main()
