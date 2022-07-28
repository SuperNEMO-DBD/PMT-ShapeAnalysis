import sys

sys.path.insert(1, '../..')
from pmt_he_study.models import *


def plot_corrected(data):
    t_popt, t_pcov, t_perr = adjust_pars()
    nums = [data[0]["he_ap_nums"], data[1]["he_ap_nums"]]
    nums_err = [data[0]["he_ap_nums_err"], data[1]["he_ap_nums_err"]]

    print(nums_err[0]/nums[0])

    err_p0 = t_perr[0] / t_popt[0]
    err_p1 = t_perr[1] / t_popt[1]
    err_p2 = t_perr[2]
    A = [t_popt[0] * (data[i]["he_ap_nums"] ** 2) for i in range(2)]
    err_A = [A[i] * np.sqrt((err_p0) ** 2 + 2 * (data[i]["he_ap_nums_err"]/data[i]["he_ap_nums"]) ** 2) for i in range(2)]
    B = [t_popt[1] * data[i]["he_ap_nums"] for i in range(2)]
    err_B = [B[i] * np.sqrt((err_p1) ** 2 + 2 * (data[i]["he_ap_nums_err"]/data[i]["he_ap_nums"]) ** 2) for i in range(2)]
    C = t_popt[2]
    err_C = err_p2
    he_nums_corr = [A[i] + B[i] + C for i in range(2)]
    he_nums_corr_err = [np.sqrt((err_A[i]) ** 2 + (err_B[i]) ** 2 + (err_C) ** 2) for i in range(2)]

    plt.figure(figsize=figsize)
    plt.errorbar(data[0]["dates"], nums[0], yerr=nums_err[0], fmt="o", label='Raw Data', markersize=1,
                 capsize=cap_size, linewidth=line_width, capthick=cap_thick)
    plt.errorbar(data[0]["dates"], he_nums_corr[0], yerr=he_nums_corr_err[0], fmt="^", label='Corrected',
                 markersize=1, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
    plt.xlabel("Days from onset of 1% Helium")
    plt.ylabel("Number")
    # plt.title("Average After-Pulse Number Correction")
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/PMT_Project/data_correction.pdf")


def read_gain_file(filename: str):
    # According to gain_plotter.py
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
    data = {
        1000: {
            0: {
                "week": np.array([]),
                "mev_charge": np.array([]),
                "av_charge": np.array([]),
                "gain": np.array([]),
                "res": np.array([]),
                "mev_charge_err": np.array([]),
                "av_charge_err": np.array([]),
                "gain_err": np.array([]),
                "res_err": np.array([]),
            },
            1: {
                "week": np.array([]),
                "mev_charge": np.array([]),
                "av_charge": np.array([]),
                "gain": np.array([]),
                "res": np.array([]),
                "mev_charge_err": np.array([]),
                "av_charge_err": np.array([]),
                "gain_err": np.array([]),
                "res_err": np.array([]),
            },
        },
        1400: {
            0: {
                "week": np.array([]),
                "mev_charge": np.array([]),
                "av_charge": np.array([]),
                "gain": np.array([]),
                "res": np.array([]),
                "mev_charge_err": np.array([]),
                "av_charge_err": np.array([]),
                "gain_err": np.array([]),
                "res_err": np.array([]),
            },
            1: {
                "week": np.array([]),
                "mev_charge": np.array([]),
                "av_charge": np.array([]),
                "gain": np.array([]),
                "res": np.array([]),
                "mev_charge_err": np.array([]),
                "av_charge_err": np.array([]),
                "gain_err": np.array([]),
                "res_err": np.array([]),
            },
        }
    }

    file = open(filename, "r")
    fl = file.readlines()
    for index, line in enumerate(fl):
        if index == 0:
            continue
        ll = line.split(",")
        if len(ll) < 10:
            continue
        for i, volt in enumerate([1000, 1400]):
            for ch in [0, 1]:
                data[volt][ch]["week"]           = np.append(data[volt][ch]["week"],            int(ll[0]))
                data[volt][ch]["mev_charge"]     = np.append(data[volt][ch]["mev_charge"],      float(ll[1 + i*4 + ch*2]))
                data[volt][ch]["mev_charge_err"] = np.append(data[volt][ch]["mev_charge_err"],  float(ll[2 + i*4 + ch*2]))
                data[volt][ch]["av_charge"]      = np.append(data[volt][ch]["av_charge"],       float(ll[9 + i*4 + ch*2]))
                data[volt][ch]["av_charge_err"]  = np.append(data[volt][ch]["av_charge_err"],   float(ll[10 + i*4 + ch*2]))
                data[volt][ch]["gain"]           = np.append(data[volt][ch]["gain"],            float(ll[17 + i*4 + ch*2]))
                data[volt][ch]["gain_err"]       = np.append(data[volt][ch]["gain_err"],        float(ll[18 + i*4 + ch*2]))
                if i == 0:
                    try:
                        data[volt][ch]["res"] = np.append(data[volt][ch]["res"], float(ll[25 + i*4 + ch*2]))
                        data[volt][ch]["res_err"] = np.append(data[volt][ch]["res_err"], float(ll[26 + i*4 + ch*2]))
                    except ValueError:
                        data[volt][ch]["res"] = np.append(data[volt][ch]["res"], None)
                        data[volt][ch]["res_err"] = np.append(data[volt][ch]["res_err"], None)

    return data


def read_file(pmt_array: PMT_Array, input_filename: str):
    try:
        print(">>> Reading data from file: {}".format(input_filename))
        date_file = open(input_filename, 'r')
    except FileNotFoundError as fnf_error:
        print(fnf_error)
        raise Exception("Error opening data file {}".format(input_filename))

    filenames = np.loadtxt(input_filename, delimiter=',', dtype={
        'names': ['filename'],
        'formats': ['S100']}, unpack=True)

    data = {
        0: {
            "dates": np.array([], dtype='int'),
            "ap_nums": np.array([]),
            "ap_nums_err": np.array([]),
            "he_ap_nums": np.array([]),
            "he_ap_nums_err": np.array([])
        },
        1: {
            "dates": np.array([], dtype='int'),
            "ap_nums": np.array([]),
            "ap_nums_err": np.array([]),
            "he_ap_nums": np.array([]),
            "he_ap_nums_err": np.array([])
        }
    }

    for i_file in range(filenames.size):
        filename = filenames[i_file][0].decode("utf-8")
        date = filename.split("_")[0]
        voltage = int(filename.split("_")[1].split("A")[1])

        file = ROOT.TFile(input_filename.split("filenames.txt")[0] + filename, "READ")
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

            if he_num_hist.GetMean() < 0.01:
                pass
            else:
                data[i_om]["dates"] = np.append(data[i_om]["dates"], process_date([int(date)]))
                data[i_om]["ap_nums"] = np.append(data[i_om]["ap_nums"], num_hist.GetMean())
                data[i_om]["ap_nums_err"] = np.append(data[i_om]["ap_nums_err"], num_hist.GetMeanError())
                data[i_om]["he_ap_nums"] = np.append(data[i_om]["he_ap_nums"], he_num_hist.GetMean())
                data[i_om]["he_ap_nums_err"] = np.append(data[i_om]["he_ap_nums_err"], he_num_hist.GetMeanError())

            del num_hist
            del he_num_hist
        file.Close()
    return data


def plot_av_charge_gain(gain_data: dict, voltage: int):
    fig1 = plt.figure(figsize=(4.5, 4), facecolor='w')
    frame1 = fig1.add_axes((.15, .54, .8, .39))
    frame1.set_xticklabels([])
    val = 1
    if voltage == 1400:
        plt.ylim(100, 150)
        val = 1e6

    plt.errorbar(gain_data[voltage][0]["week"] * 7 + 3, gain_data[voltage][0]["av_charge"],
                 yerr=gain_data[voltage][0]["av_charge_err"],
                 fmt="o", label='Exposed', markersize=1, capsize=cap_size, linewidth=line_width,
                 capthick=cap_thick)
    plt.errorbar(gain_data[voltage][1]["week"] * 7 + 3, gain_data[voltage][1]["av_charge"],
                 yerr=gain_data[voltage][1]["av_charge_err"],
                 fmt="^", label='Control', markersize=1, capsize=cap_size, linewidth=line_width,
                 capthick=cap_thick)
    plt.legend(loc='best')
    plt.ylabel("Charge /pC")

    frame2 = fig1.add_axes((.15, .12, .8, 0.39))

    plt.errorbar(gain_data[voltage][0]["week"] * 7 + 3, gain_data[voltage][0]["gain"]/val,
                 yerr=gain_data[voltage][0]["gain_err"]/val,
                 fmt="o", label='Exposed', markersize=1, capsize=cap_size, linewidth=line_width,
                 capthick=cap_thick)
    plt.errorbar(gain_data[voltage][1]["week"] * 7 + 3, gain_data[voltage][1]["gain"]/val,
                 yerr=gain_data[voltage][1]["gain_err"]/val,
                 fmt="^", label='Control', markersize=marker_size, capsize=cap_size, linewidth=line_width,
                 capthick=cap_thick)
    plt.xlabel("Days from onset of 1% Helium")
    plt.ylabel(r"Gain /$10^6$")
    # ax2.set_ylim(1, 4)

    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/PMT_Project/av_charge_gain_{}V.pdf".format(voltage))


def plot_aan(data: dict):
    t_popt, t_pcov, t_perr = adjust_pars()
    try:
        start = np.where(data[0]["dates"] == 0)[0][0]
    except:
        start = np.where(data[0]["dates"] == 1)[0][0]
    mid = np.where(data[0]["dates"] == 98)[0][0]

    nums = [data[0]["ap_nums"], data[1]["ap_nums"]]
    nums_err = [data[0]["ap_nums_err"], data[1]["ap_nums_err"]]
    err_p0 = t_perr[0] / t_popt[0]
    err_p1 = t_perr[1] / t_popt[1]
    err_p2 = t_perr[2]
    A = [t_popt[0] * (data[i]["he_ap_nums"] ** 2) for i in range(2)]
    err_A = [A[i] * np.sqrt((err_p0) ** 2 + 2 * (data[i]["he_ap_nums_err"]/data[i]["he_ap_nums"]) ** 2) for i in range(2)]
    B = [t_popt[1] * data[i]["he_ap_nums"] for i in range(2)]
    err_B = [B[i] * np.sqrt((err_p1) ** 2 + 2 * (data[i]["he_ap_nums_err"]/data[i]["he_ap_nums"]) ** 2) for i in range(2)]
    C = t_popt[2]
    err_C = err_p2
    he_nums_corr = [A[i] + B[i] + C for i in range(2)]
    he_nums_corr_err = [np.sqrt((err_A[i]) ** 2 + (err_B[i]) ** 2 + (err_C) ** 2) for i in range(2)]

    x_0, y_0, y_err_0 = data[0]["dates"][:start + 1], nums[0][:start + 1], nums_err[0][:start + 1]
    x_1, y_1, y_err_1 = data[0]["dates"][start + 1:mid + 1], nums[0][start + 1:mid + 1], nums_err[0][start + 1:mid + 1]
    x_2, y_2, y_err_2 = data[0]["dates"][mid + 1:], nums[0][mid + 1:], nums_err[0][mid + 1:]
    x_3, y_3, y_err_3 = data[1]["dates"], nums[1], nums_err[1]

    plt.figure(figsize=figsize, facecolor='w')
    plt.errorbar(x_0, y_0, zorder=0, yerr=y_err_0, fmt="C0s", label="Atmospheric He",
                 markersize=1, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
    plt.errorbar(x_1, y_1, zorder=0, yerr=y_err_1, fmt="C1s", label="1% He",
                 markersize=1, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
    plt.errorbar(x_2, y_2, zorder=0, yerr=y_err_2, fmt="C2s", label="10% He",
                 markersize=1, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
    plt.errorbar(x_3, y_3, zorder=0, yerr=y_err_3, fmt="C3o", label="Control",
                 markersize=1, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
    plt.ylabel("Number")
    plt.title("Average After-pulse Number")
    plt.xlim(-30, 420)
    plt.xlabel("Days from 1% Helium Onset")
    plt.legend(loc='upper left')
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/PMT_Project/aan_vs_time.pdf")

    x_0, y_0, y_err_0 = data[0]["dates"][:start + 1], he_nums_corr[0][:start + 1], he_nums_corr_err[0][:start + 1]
    x_1, y_1, y_err_1 = data[0]["dates"][start + 1:mid + 1], he_nums_corr[0][start + 1:mid + 1], he_nums_corr_err[0][start + 1:mid + 1]
    x_2, y_2, y_err_2 = data[0]["dates"][mid + 1:], he_nums_corr[0][mid + 1:], he_nums_corr_err[0][mid + 1:]
    x_3, y_3, y_err_3 = data[1]["dates"], he_nums_corr[1], he_nums_corr_err[1]

    plt.figure(figsize=figsize, facecolor='w')
    plt.errorbar(x_0, y_0, zorder=0, yerr=y_err_0, fmt="C0s", label="Atmospheric He",
                 markersize=1, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
    plt.errorbar(x_1, y_1, zorder=0, yerr=y_err_1, fmt="C1s", label="1% He",
                 markersize=1, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
    plt.errorbar(x_2, y_2, zorder=0, yerr=y_err_2, fmt="C2s", label="10% He",
                 markersize=1, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
    plt.errorbar(x_3, y_3, zorder=0, yerr=y_err_3, fmt="C3o", label="Control",
                 markersize=1, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
    plt.errorbar(data[0]["dates"], data[0]["he_ap_nums"], zorder=0, yerr=data[0]["he_ap_nums_err"], fmt="C4o",
                 label="Exposed raw",
                 markersize=1, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
    plt.errorbar(data[1]["dates"], data[1]["he_ap_nums"], zorder=0, yerr=data[1]["he_ap_nums_err"], fmt="C5o",
                 label="Control raw",
                 markersize=1, capsize=cap_size, linewidth=line_width, capthick=cap_thick)

    plt.ylabel("Number")
    plt.title(r"He$^+$ Average After-pulse Number")
    plt.xlim(-30, 420)
    plt.xlabel("Days from 1% Helium Onset")
    plt.legend(loc='upper left')
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/PMT_Project/he_aan_vs_time.pdf")


def plot_pp(model, data: dict, gain_data: dict):
    t_popt, t_pcov, t_perr = adjust_pars()
    try:
        start = np.where(data[0]["dates"] == 0)[0][0]
    except:
        start = np.where(data[0]["dates"] == 1)[0][0]
    mid = np.where(data[0]["dates"] == 98)[0][0]

    nums = [data[0]["ap_nums"], data[1]["ap_nums"]]
    # nums_err = [data[0]["ap_nums"] * 0.01, data[1]["ap_nums"] * 0.01]
    nums_err = [data[0]["ap_nums_err"], data[1]["ap_nums_err"]]

    err_p0 = t_perr[0] / t_popt[0]
    err_p1 = t_perr[1] / t_popt[1]
    err_p2 = t_perr[2]
    A = [t_popt[0] * (data[i]["he_ap_nums"] ** 2) for i in range(2)]
    err_A = [A[i]*np.sqrt((err_p0)**2 + 2*(data[i]["he_ap_nums_err"]/data[i]["he_ap_nums"])**2) for i in range(2)]
    B = [t_popt[1] * data[i]["he_ap_nums"] for i in range(2)]
    err_B = [B[i]*np.sqrt((err_p1)**2 + 2*(data[i]["he_ap_nums_err"]/data[i]["he_ap_nums"])**2) for i in range(2)]
    C = t_popt[2]
    err_C = err_p2
    he_nums_corr = [A[i] + B[i] + C for i in range(2)]
    he_nums_corr_err = [np.sqrt((err_A[i])**2 + (err_B[i])**2 + (err_C)**2) for i in range(2)]
    pp, pp_err = partial_pressure(data, gain_data, nums, nums_err)
    he_pp, he_pp_err = partial_pressure(data, gain_data, he_nums_corr, he_nums_corr_err)

    x = data[0]["dates"][mid + 1:] - data[0]["dates"][mid + 1:][0]
    y = pp[0][mid + 1:]
    yerr = pp_err[0][mid + 1:]
    pars, errs, chi, ndof = root_fit(x, y, yerr, model)

    x_model = np.array(data[0]["dates"][mid + 1:])
    y_model = model.func(x, pars)

    x_0, y_0, y_err_0 = data[0]["dates"][:start + 1],        pp[0][:start + 1],         pp_err[0][:start + 1]
    x_1, y_1, y_err_1 = data[0]["dates"][start + 1:mid + 1], pp[0][start + 1:mid + 1],  pp_err[0][start + 1:mid + 1]
    x_2, y_2, y_err_2 = data[0]["dates"][mid + 1:],          pp[0][mid + 1:],           pp_err[0][mid + 1:]
    x_3, y_3, y_err_3 = data[1]["dates"],                    pp[1],                     pp_err[1]

    k = pars[1]
    # extrapolate(model, pars, 3.13, "")
    # extrapolate(model, pars, 0.13, "")

    fig1 = plt.figure(figsize=(5, 4), facecolor='w')
    frame1 = fig1.add_axes((.15, .32, .8, .6))
    frame1.set_xticklabels([])
    plt.errorbar(x_0, y_0 - k, zorder=0, yerr=y_err_0, fmt="C0s", label="Atmospheric He",
                 markersize=1, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
    plt.errorbar(x_1, y_1 - k, zorder=0, yerr=y_err_1, fmt="C1s", label="1% He",
                 markersize=1, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
    plt.errorbar(x_2, y_2 - k, zorder=0, yerr=y_err_2, fmt="C2s", label="10% He",
                 markersize=1, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
    plt.errorbar(x_3, y_3 - np.average(y_3[mid + 1:]), zorder=0, yerr=y_err_3, fmt="C3o", label="Control",
                 markersize=1, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
    plt.plot(x_model, y_model - k, 'k-', label='Model', zorder=10, linewidth=line_width)

    handles, labels = plt.gca().get_legend_handles_labels()

    n_sf = [0, 0, 0]
    for i in range(3):
        while errs[i]/pars[i] < 1:
            if errs[i]/pars[i] * pow(10, n_sf[i]) >= 1:
                break
            else:
                n_sf[i] += 1
    strings = [r'$P_0 =$ {:.' + str(n_sf[0]) + 'e} ± {:.0e}', r'$P_1 =$ {:.' + str(n_sf[1]) + 'e} ± {:.0e}',
               '$L =$ {:.0f} ± {:.0f} days']
    patch = patches.Patch(color='white', label=strings[0].format(pars[0], errs[0]))
    patch_1 = patches.Patch(color='white', label=strings[1].format(pars[1], errs[1]))
    patch_2 = patches.Patch(color='white',
                            label=strings[2].format(pars[2] / (3600 * 24), errs[2] / (3600 * 24)))
    patch_3 = patches.Patch(color='white', label=r'$\chi^2$/N$_{DoF}$' + ' = {:.2f}/{}'.format(chi*ndof, ndof))
    handles.extend([patch, patch_1, patch_2, patch_3])

    plt.ylabel("Pressure /Pa")
    # plt.title("Reconstructed Helium Internal Partial Pressure")
    plt.xlim(-30, 420)
    plt.legend(handles=handles, loc='upper left')

    frame2 = fig1.add_axes((.15, .1, .8, .2))
    plt.xlabel("Days from 1% Helium Onset")
    plt.axhline(0, ls='--', color='black')
    plt.ylabel("(model-data)/model")
    plt.xlim(-30, 420)
    plt.errorbar(x_model, (y_model - y_2) / y_model, yerr=y_err_2 / y_model, fmt="k.",
                 markersize=1, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
    plt.tight_layout()

    plt.savefig("/Users/williamquinn/Desktop/PMT_Project/" + model.name + "_pi_vs_time.pdf")

    x = data[0]["dates"][mid + 1:] - data[0]["dates"][mid + 1:][0]
    y = he_pp[0][mid + 1:]
    yerr = he_pp_err[0][mid + 1:]
    pars, errs, chi, ndof = root_fit(x, y, yerr, model)

    x_model = np.array(data[0]["dates"][mid + 1:])
    y_model = model.func(x, pars)

    x_0, y_0, y_err_0 = data[0]["dates"][:start + 1], he_pp[0][:start + 1], he_pp_err[0][:start + 1]
    x_1, y_1, y_err_1 = data[0]["dates"][start + 1:mid + 1], he_pp[0][start + 1:mid + 1], he_pp_err[0][start + 1:mid + 1]
    x_2, y_2, y_err_2 = data[0]["dates"][mid + 1:], he_pp[0][mid + 1:], he_pp_err[0][mid + 1:]
    x_3, y_3, y_err_3 = data[1]["dates"], he_pp[1], he_pp_err[1]

    k = pars[1]
    # extrapolate(model, pars, 3.13, "he")
    # extrapolate(model, pars, 0.13, "he")

    fig1 = plt.figure(figsize=(5, 4), facecolor='w')
    frame1 = fig1.add_axes((.15, .32, .8, .6))
    frame1.set_xticklabels([])
    plt.errorbar(x_0, y_0 - k, zorder=0, yerr=y_err_0, fmt="C0s", label="Atmospheric He",
                 markersize=1, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
    plt.errorbar(x_1, y_1 - k, zorder=0, yerr=y_err_1, fmt="C1s", label="1% He",
                 markersize=1, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
    plt.errorbar(x_2, y_2 - k, zorder=0, yerr=y_err_2, fmt="C2s", label="10% He",
                 markersize=1, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
    plt.errorbar(x_3, y_3 - np.average(y_3[mid + 1:]), zorder=0, yerr=y_err_3, fmt="C3o", label="Control",
                 markersize=1, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
    plt.plot(x_model, y_model - k, 'k-', label='Model', zorder=10, linewidth=line_width)

    handles, labels = plt.gca().get_legend_handles_labels()
    n_sf = [0, 0, 0]
    for i in range(3):
        while 1:
            if (errs[i] / pars[i]) * pow(10, n_sf[i]) >= 1:
                break
            else:
                n_sf[i] += 1
    strings = [r'$P_0 =$ {:.' + str(n_sf[0]) + 'e} ± {:.0e}', r'$P_1 =$ {:.' + str(n_sf[1]) + 'e} ± {:.0e}',
               '$L =$ {:.0f} ± {:.0f} days']
    patch = patches.Patch(color='white', label=strings[0].format(pars[0], errs[0]))
    patch_1 = patches.Patch(color='white', label=strings[1].format(pars[1], errs[1]))
    patch_2 = patches.Patch(color='white',
                            label=strings[2].format(pars[2] / (3600 * 24), errs[2] / (3600 * 24)))
    patch_3 = patches.Patch(color='white', label=r'$\chi^2$/N$_{DoF}$' + ' = {:.2f}/{}'.format(chi*ndof, ndof))
    handles.extend([patch, patch_1, patch_2, patch_3])

    plt.ylabel("Pressure /Pa")
    # plt.title("Reconstructed Helium Internal Partial Pressure")
    plt.xlim(-30, 420)
    plt.legend(handles=handles, loc='upper left')

    frame2 = fig1.add_axes((.15, .1, .8, .2))
    plt.xlabel("Days from 1% Helium Onset")
    plt.axhline(0, ls='--', color='black')
    plt.ylabel("(model-data)/model")
    plt.xlim(-30, 420)
    plt.errorbar(x_model, (y_model - y_2) / y_model, yerr=y_err_2 / y_model, fmt="k.",
                 markersize=1, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
    plt.tight_layout()

    plt.savefig("/Users/williamquinn/Desktop/PMT_Project/" + model.name + "_he_pi_vs_time.pdf")


def main():
    gain_data = read_gain_file("/Users/williamquinn/Desktop/PMT_Project/week_values.csv")

    topology = [2, 1]
    pmt_array = PMT_Array(topology, "summary")
    pmt_array.set_pmt_id("GAO612", 1)
    pmt_array.set_pmt_id("GAO607", 0)

    data = read_file(pmt_array, "/Users/williamquinn/Desktop/PMT_Project/data/1400V/filenames.txt")
    # plot_aan(data)
    # plot_corrected(data)
    # plot_av_charge_gain(gain_data, 1400)

    model = Model()
    plot_pp(model, data, gain_data)
    model = Model_0()
    plot_pp(model, data, gain_data)


if __name__ == "__main__":
    main()
