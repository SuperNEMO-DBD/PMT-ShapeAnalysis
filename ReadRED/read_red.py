import ROOT
import numpy as np
import matplotlib.pyplot as plt
import sndisplay as sn
from format_plot import *
from scipy.optimize import curve_fit
from random import randint


c_tdc2sec = 0.625E-08
t_tdc2sec = 1.25E-08
tdc2ns = 0.390625
adc2mv = 0.610352
MAX_GG_TIMES = 10


def cell_id(cellnum):
    cell_side = cellnum // (9 * 113)
    cell_row = cellnum % (9 * 113) // 9
    cell_layer = cellnum % (9 * 113) % 9

    return cell_side, cell_row, cell_layer


def gaussian(x, a, mu, sig):
    return a*np.exp(-0.5*((x - mu)/sig)**2)


def gaussian_B(x, a, mu, sig, b):
    return a * np.exp(-0.5 * ((x - mu) / sig) ** 2) + b


def parse_arguments():
    import argparse
    parser = argparse.ArgumentParser(description="Input file names")
    parser.add_argument('-i', required=True, type=str, help='Input data file')
    parser.add_argument('-d', required=False, type=str, help='Output Directory')
    parser.add_argument('-time', required=False, type=str, help='run time')
    args = parser.parse_args()
    return args


def plot_bad_map(anode_events, top_cathode_events, bottom_cathode_events, output_directory):
    sntracker = sn.tracker("bad_map", with_palette=True)
    sntracker.draw_cellid_label()
    sntracker.draw_content_label('{}')
    na = 0
    nb = 0
    nt = 0
    nboth = 0

    for i_cell in range(2034):
        bad_anode = False
        bad_t_c = False
        bad_b_c = False
        if anode_events[i_cell] == 0 or anode_events[i_cell] < 20:
            bad_anode = True
        if top_cathode_events[i_cell] == 0 or top_cathode_events[i_cell] < 10:
            bad_t_c = True
        if bottom_cathode_events[i_cell] == 0 or bottom_cathode_events[i_cell] < 10:
            bad_b_c = True

        if bad_anode:
            sntracker.setcontent(i_cell, 1)
            na += 1
        elif bad_t_c and bad_b_c:
            sntracker.setcontent(i_cell, 2)
            nboth += 1
        elif bad_t_c:
            sntracker.setcontent(i_cell, 3)
            nt += 1
        elif bad_b_c:
            sntracker.setcontent(i_cell, 4)
            nb += 1

    print("a:", na, "both:", nboth, "nb:", nb, "nt:", nt)
    # sntracker.setrange(0, 35)

    sntracker.draw()
    sntracker.save(output_directory)


def plot_cathode_events(tree: ROOT.TTree, cut: bool, run, output_directory, run_time):
    top_cathode_events = [0 for i in range(2034)]
    bottom_cathode_events = [0 for i in range(2034)]

    if cut:
        top_name = run + "_top_cathode_events"
        bot_name = run + "_bottom_cathode_events"
    else:
        top_name = run + "_top_cathode_events_cut"
        bot_name = run + "_bottom_cathode_events_cut"

    for event in tree:
        if cut:
            if len(event.tracker_time_anode) > 50:
                continue
            if len(event.calo_time) == 0:
                continue
            if len(event.tracker_time_anode) != len(event.tracker_time_top_cathode):
                continue
            times = []
            for j in range(len(event.calo_time)):
                if bool(event.calo_high_t[j]):
                    times.append(event.calo_time[j])
            if len(times) == 0:
                continue

        # print(len(event.tracker_cell_num), len(event.tracker_time_top_cathode))

        for index, i_cell in enumerate(event.tracker_cell_num):
            if event.tracker_time_top_cathode[index] == is_none:
                pass
            else:
                top_cathode_events[event.tracker_cell_num[index]] += 1/run_time
            if event.tracker_time_bottom_cathode[index] == is_none:
                pass
            else:
                bottom_cathode_events[event.tracker_cell_num[index]] += 1/run_time
    plot_event_map(top_name, top_cathode_events, output_directory)
    plot_event_map(bot_name, bottom_cathode_events, output_directory)


def plot_anode_events(tree: ROOT.TTree, cut: bool, run, output_directory, run_time):
    anode_events = [0 for i in range(2034)]

    if cut:
        name = run + "_anode_events"
    else:
        name = run + "_anode_events_cut"

    for event in tree:
        if cut:
            if len(event.tracker_time_anode) > 50:
                continue
            if len(event.calo_time) == 0:
                continue
            if len(event.tracker_time_anode) != len(event.tracker_time_top_cathode):
                continue
            times = []
            for j in range(len(event.calo_time)):
                if bool(event.calo_high_t[j]):
                    times.append(event.calo_time[j])
            if len(times) == 0:
                continue

        # print(len(event.tracker_cell_num), len(event.tracker_time_top_cathode))

        for index, i_cell in enumerate(event.tracker_cell_num):
            if event.tracker_time_anode[index] == is_none:
                pass
            else:
                anode_events[event.tracker_cell_num[index]] += 1 / run_time
    plot_event_map(name, anode_events, output_directory)


def plot_drift_times_per_cell(drifts, name, run_time, output_directory):
    means = []
    stds = []
    chis = []
    lower, upper = 0, 100
    n_bins = 100
    for i_cell in range(len(drifts)):
        if len(drifts[i_cell]) > 0:
            mean = np.average(drifts[i_cell])
            std = np.std(drifts[i_cell])
            freq, bin_edges = np.histogram(drifts[i_cell], n_bins, range=(lower, upper))
            width = bin_edges[-1] - bin_edges[-2]
            bin_centres = bin_edges[:-1] + width / 2

            fit = ROOT.TF1("fit", "[0]*TMath::Gaus(x, [1], [2])", lower, upper)
            hist = ROOT.TH1D("", "", n_bins, lower, upper)
            for i in range(len(drifts[i_cell])):
                hist.Fill(drifts[i_cell][i])
            hist.Scale(1/run_time)
            fit.SetParLimits(0, 0, (np.max(freq) + np.max(freq)*0.1)/run_time)
            fit.SetParLimits(1, bin_centres[np.argmax(freq)] - 10, bin_centres[np.argmax(freq)] + 10)
            fit.SetParLimits(2, 0, 50)
            fit.SetParameters(np.max(freq), bin_centres[np.argmax(freq)], 1)
            hist.Fit("fit", "0Q")
            A = fit.GetParameter(0)
            mu = fit.GetParameter(1)
            sig = fit.GetParameter(2)
            x = np.linspace(lower, upper, 500)

            fig = plt.figure(figsize=figsize, facecolor='white')
            plt.bar(bin_centres, freq/run_time, width=width)
            plt.plot(x, gaussian(x, A, mu, sig), 'C3', label='model')
            plt.axvline(mu, ls='--', color='C1', label='Mean')
            plt.axvline(mu - sig, ls='--', color='C2')
            plt.axvline(mu + sig, ls='--', color='C2')
            plt.title('Cell: {}.{}.{} Mean: {:.2f}'.format(*cell_id(i_cell), mu))
            plt.xlabel("Drift Time /µs")
            plt.ylabel("Counts per Minute")
            plt.legend(loc='best')
            plt.xlim(lower, upper)
            plt.tight_layout()
            plt.savefig(output_directory + name + '_cell_{}-{}-{}'.format(*cell_id(i_cell)) + '.pdf')
            plt.close()

            if fit.GetNDF() == 0:
                chi = None
            else:
                chi = fit.GetChisquare() / fit.GetNDF()

            means.append(mu)
            stds.append(sig)
            chis.append(chi)

            del hist
            del fit

        else:
            means.append(None)
            stds.append(None)
            chis.append(None)

    sntracker = sn.tracker(name + "_av", True)
    sntracker.draw_cellid_label()
    sntracker.draw_content_label('{:.2f}')

    for i_cell in range(len(means)):
        if means[i_cell] is None:
            continue
        sntracker.setcontent(i_cell, means[i_cell])

    sntracker.setrange(30, 80)

    sntracker.draw()
    sntracker.save(output_directory)

    del sntracker

    sntracker = sn.tracker(name + "_std", True)
    sntracker.draw_cellid_label()
    sntracker.draw_content_label('{:.2f}')

    for i_cell in range(len(stds)):
        if stds[i_cell] is None:
            continue
        sntracker.setcontent(i_cell, stds[i_cell])

    sntracker.draw()
    sntracker.save(output_directory)

    del sntracker

    sntracker = sn.tracker(name + "_chi2", True)
    sntracker.draw_cellid_label()
    sntracker.draw_content_label('{:.2f}')

    for i_cell in range(len(stds)):
        if chis[i_cell] is None:
            continue
        sntracker.setcontent(i_cell, chis[i_cell])

    sntracker.draw()
    sntracker.save(output_directory)

    del sntracker


def plot_drift_times_all(event_tree: ROOT.TTree, run_time, output_directory, run_num):
    anode_times = []
    first_cathode = []
    second_cathode = []
    top_cathode_times = []
    bot_cathode_times = []
    drifts = [[] for i in range(2034)]
    drifts_0 = [[] for i in range(2034)]
    n_event = 0
    for event in event_tree:
        # Cut out an event that has more that 50 tracker hits
        if len(event.tracker_time_anode) > 50:
            continue
        '''# Cut out events that have no calorimeter hits
        # TODO: treat this case seperatly
        if len(event.calo_time) == 0:
            continue'''
        calo_time = None
        ht = False
        if len(event.calo_time) != 0:
            for j in range(len(event.calo_time)):
                if bool(event.calo_high_t[j]):
                    ht = True
                    temp_calo_time = event.calo_time[j]
                    if calo_time is not None and temp_calo_time < calo_time:
                        calo_time = temp_calo_time
                    elif calo_time is None:
                        calo_time = temp_calo_time

        for i in range(len(event.tracker_time_anode)):
            if event.tracker_time_top_cathode[i] != is_none and event.tracker_time_bottom_cathode[i] != is_none and event.tracker_time_anode[i] != is_none:
                time = event.tracker_time_anode[i]
                if calo_time is not None:
                    anode_times.append((event.tracker_time_anode[i] - calo_time) * 1e6)
                top_cathode_times.append((event.tracker_time_top_cathode[i] - time) * 1e6)
                bot_cathode_times.append((event.tracker_time_bottom_cathode[i] - time) * 1e6)
                first_cathode.append((event.tracker_time_anode_first_lt[i] - time) * 1e6)
                second_cathode.append((event.tracker_time_anode_second_lt[i] - time) * 1e6)

                drift = ((event.tracker_time_top_cathode[i] - time) + (event.tracker_time_bottom_cathode[i] - time)) * 1e6
                drift_0 = ((event.tracker_time_anode_first_lt[i] - time) + (event.tracker_time_anode_second_lt[i] - time)) * 1e6
                drifts[event.tracker_cell_num[i]].append(drift)
                drifts_0[event.tracker_cell_num[i]].append(drift_0)

        n_event += 1

    lower, upper = -50, 200
    fig = plt.figure(figsize=figsize, facecolor='white')
    freq, bin_edges = np.histogram(np.array(anode_times), 250, range=(lower, upper))
    width = bin_edges[-1] - bin_edges[-2]
    bin_centres = bin_edges[:-1] + width / 2
    plt.plot(bin_centres, freq/run_time, "-", label='anode')

    freq, bin_edges = np.histogram(np.array(top_cathode_times), 250, range=(lower, upper))
    width = bin_edges[-1] - bin_edges[-2]
    bin_centres = bin_edges[:-1] + width / 2
    plt.plot(bin_centres, freq/run_time, "-", label='top cathode')

    freq, bin_edges = np.histogram(np.array(bot_cathode_times), 250, range=(lower, upper))
    width = bin_edges[-1] - bin_edges[-2]
    bin_centres = bin_edges[:-1] + width / 2
    plt.plot(bin_centres, freq/run_time, "-", label='bottom cathode')

    freq, bin_edges = np.histogram(np.array(first_cathode), 250, range=(lower, upper))
    width = bin_edges[-1] - bin_edges[-2]
    bin_centres = bin_edges[:-1] + width / 2
    plt.plot(bin_centres, freq/run_time, "-", label='first cathode')

    freq, bin_edges = np.histogram(np.array(second_cathode), 250, range=(lower, upper))
    width = bin_edges[-1] - bin_edges[-2]
    bin_centres = bin_edges[:-1] + width / 2
    plt.plot(bin_centres, freq/run_time, "-", label='second cathode')

    plt.xlabel("Drift Time /µs")
    plt.ylabel("Counts per minute")
    plt.title("Run {} All Cells".format(run_num))
    plt.yscale('log')
    plt.legend(loc='best')
    plt.xlim(lower, upper)
    plt.tight_layout()
    plt.savefig(output_directory + "/plots/run_" + str(run_num) + "/" + str(run_num) + "_all_times.pdf")
    plt.close()

    all_drifts = []
    for i_cell in range(len(drifts)):
        for i_val in range(len(drifts[i_cell])):
            all_drifts.append(drifts[i_cell][i_val])

    fig = plt.figure(figsize=figsize, facecolor='white')
    freq, bin_edges = np.histogram(np.array(all_drifts), 100, range=(0, 100))
    width = bin_edges[-1] - bin_edges[-2]
    bin_centres = bin_edges[:-1] + width / 2
    plt.bar(bin_centres, freq / run_time, width=width)

    plt.xlabel("Drift Time /µs")
    plt.ylabel("Counts per minute")
    plt.title("Run {} All Cells".format(run_num))
    plt.xlim(0, 100)
    plt.tight_layout()
    plt.savefig(output_directory + "/plots/run_" + str(run_num) + "/" + str(run_num) + "_all_drift_times.pdf")
    plt.close()

    plot_drift_times_per_cell(drifts, str(run_num) + "_cell_drift_time", run_time, output_directory + "/plots/run_" + str(run_num) + "/")
    # plot_drift_times_per_cell(drifts_0, str(run_num) + "_cell_drift_time_fa", run_time, output_directory + "/plots/run_" + str(run_num) + "/")


def plot_full_event_map(directory: str, cut: bool, output_directory):
    files = ["red_617_output.root", "red_635_output.root", "red_638_output.root", "red_645_output.root",
             "red_652_output.root", "red_660_output.root", "red_668_output.root", "red_608_output.root"]
    run_times = [10, 3.666667, 30, 15, 10, 15, 10, 60]
    top_cathode_events = [0 for i in range(2034)]
    bottom_cathode_events = [0 for i in range(2034)]
    anode_events = [0 for i in range(2034)]

    bad_map = [0 for i in range(2034)]

    if cut:
        top_name = "top_cathode_events_cut"
        bot_name = "bottom_cathode_events_cut"
        anode_name = "anode_events_cut"
    else:
        top_name = "top_cathode_events"
        bot_name = "bottom_cathode_events"
        anode_name = "anode_events"

    for i_file, file in enumerate(files):
        root_file = ROOT.TFile(directory + "/" + file, "READ")
        event_tree = root_file.event_tree
        run_time = run_times[i_file]

        temp_t_c = [0 for i in range(2034)]
        temp_b_c = [0 for i in range(2034)]
        temp_a = [0 for i in range(2034)]

        for event in event_tree:
            if cut:
                if len(event.tracker_time_anode) > 50:
                    continue
                if len(event.calo_time) == 0:
                    continue
                ht = False
                for j in range(len(event.calo_time)):
                    if bool(event.calo_high_t[j]):
                        ht = True
                if not ht:
                    continue

            # print(len(event.tracker_cell_num), len(event.tracker_time_top_cathode))

            for index, i_cell in enumerate(event.tracker_cell_num):
                if event.tracker_time_top_cathode[index] == is_none:
                    pass
                else:
                    temp_t_c[event.tracker_cell_num[index]] += 1
                if event.tracker_time_bottom_cathode[index] == is_none:
                    pass
                else:
                    temp_b_c[event.tracker_cell_num[index]] += 1
                if event.tracker_time_anode[index] == is_none:
                    pass
                else:
                    temp_a[event.tracker_cell_num[index]] += 1

        a_scale = np.max(temp_a)
        b_c_scale = np.max(temp_b_c)
        t_c_scale = np.max(temp_t_c)

        for i_cell in range(len(temp_a)):
            if temp_a[i_cell] == 0:
                pass
            else:
                anode_events[i_cell] += (temp_a[i_cell] / a_scale) * 100
            if temp_b_c[i_cell] == 0:
                pass
            else:
                bottom_cathode_events[i_cell] += (temp_b_c[i_cell] / b_c_scale) * 100
            if temp_t_c[i_cell] == 0:
                pass
            else:
                top_cathode_events[i_cell] += (temp_t_c[i_cell]/t_c_scale) * 100

    plot_event_map(top_name, top_cathode_events, output_directory)
    plot_event_map(bot_name, bottom_cathode_events, output_directory)
    plot_event_map(anode_name, anode_events, output_directory)
    if cut:
        plot_bad_map(anode_events, top_cathode_events, bottom_cathode_events, output_directory)


def plot_full_drift_times(directory: str, name: str, output_directory):
    files = ["red_617_output.root", "red_635_output.root", "red_638_output.root", "red_645_output.root",
             "red_652_output.root", "red_660_output.root", "red_668_output.root", "red_608_output.root"]
    drifts = [[] for i in range(2034)]
    lower, upper = -50, 200
    n_bins = 50

    for i_file, file in enumerate(files):
        root_file = ROOT.TFile(directory + "/" + file, "READ")
        event_tree = root_file.event_tree
        tot_drift_0 = []
        tot_drift_1 = []
        for event in event_tree:
            if len(event.tracker_time_anode) > 50:
                continue
            '''if len(event.calo_time) == 0:
                continue
            ht = False
            for j in range(len(event.calo_time)):
                if bool(event.calo_high_t[j]):
                    ht = True
            if not ht:
                continue'''

            for i in range(len(event.tracker_time_anode)):
                if event.tracker_time_top_cathode[i] != is_none and event.tracker_time_bottom_cathode[i] != is_none and event.tracker_time_anode[i] != is_none:
                    time = event.tracker_time_anode[i]
                    drift = ((event.tracker_time_top_cathode[i] - time) + (event.tracker_time_bottom_cathode[i] - time))*1e6
                    drifts[event.tracker_cell_num[i]].append(drift)

                    side, row, layer = cell_id(event.tracker_cell_num[i])
                    if layer in [0, 8] or row in [0, 13, 14, 27, 28, 41, 42, 55, 56, 68, 71, 84, 85, 98, 99, 112]:
                        tot_drift_0.append(drift)
                    else:
                        tot_drift_1.append(drift)

        '''fig = plt.figure(figsize=figsize, facecolor='white')
        freq, bin_edges = np.histogram(tot_drift_0, 100, range=(0, 100))
        width = bin_edges[-1] - bin_edges[-2]
        bin_centres = bin_edges[:-1] + width / 2
        plt.bar(bin_centres, freq, width=width, alpha=0.5, label='Edges')
        # popt_0, pcov_0 = curve_fit(f=gaussian_B, xdata=bin_centres, ydata=freq, sigma=np.sqrt(freq))
        # plt.plot(np.linspace(0, 100, 500), gaussian_B(np.linspace(0, 100, 500), *popt_0), 'k-')
        freq, bin_edges = np.histogram(tot_drift_1, 100, range=(0, 100))
        width = bin_edges[-1] - bin_edges[-2]
        bin_centres = bin_edges[:-1] + width / 2
        plt.bar(bin_centres, freq, width=width, alpha=0.5, label='Centres')
        # popt_1, pcov_1 = curve_fit(f=gaussian_B, xdata=bin_centres, ydata=freq, sigma=np.sqrt(freq))
        # plt.plot(np.linspace(0, 100, 500), gaussian_B(np.linspace(0, 100, 500), *popt_1), 'k-')
        plt.xlabel("Plasma Propagation Time /µs")
        plt.ylabel("Counts")
        plt.legend(loc='best')
        plt.xlim(0, 100)
        plt.tight_layout()
        plt.savefig(output_directory + "/run_{}_all_cell_drift_times.pdf".format(file.split("_")[1]))
        plt.close()'''

    means = []
    stds = []
    chis = []
    for i_cell in range(len(drifts)):
        if len(drifts[i_cell]) > 0:
            mean = np.average(drifts[i_cell])
            std = np.std(drifts[i_cell])
            freq, bin_edges = np.histogram(drifts[i_cell], n_bins, range=(lower, upper))
            width = bin_edges[-1] - bin_edges[-2]
            bin_centres = bin_edges[:-1] + width / 2

            fit = ROOT.TF1("fit", "[0]*TMath::Gaus(x, [1], [2])", lower, upper)
            hist = ROOT.TH1D("", "", n_bins, lower, upper)
            for i in range(len(drifts[i_cell])):
                hist.Fill(drifts[i_cell][i])
            fit.SetParLimits(0, 0, (np.max(freq) + np.max(freq)*0.1))
            if bin_centres[np.argmax(freq)] - 10 < 0 or bin_centres[np.argmax(freq)] + 10 > 100:
                fit.SetParLimits(1, 0, 100)
            else:
                fit.SetParLimits(1, bin_centres[np.argmax(freq)] - 10, bin_centres[np.argmax(freq)] + 10)
            fit.SetParLimits(2, 0, 50)
            fit.SetParameters(np.max(freq), bin_centres[np.argmax(freq)], 1)
            hist.Fit("fit", "0Q")
            A = fit.GetParameter(0)
            mu = fit.GetParameter(1)
            sig = fit.GetParameter(2)
            x = np.linspace(lower, upper, 500)

            if fit.GetNDF() == 0:
                chi = None
            else:
                chi = fit.GetChisquare() / fit.GetNDF()

            means.append(mu)
            stds.append(sig)
            chis.append(chi)

            del hist
            del fit

        else:
            means.append(None)
            stds.append(None)
            chis.append(None)

    sntracker = sn.tracker(name + "_av", True)
    sntracker.draw_cellid_label()
    sntracker.draw_content_label('{:.2f}')

    for i_cell in range(len(means)):
        if means[i_cell] is None:
            continue
        sntracker.setcontent(i_cell, means[i_cell])

    sntracker.setrange(30, 80)

    sntracker.draw()
    sntracker.save(output_directory)

    sntracker = sn.tracker(name + "_std", True)
    sntracker.draw_cellid_label()
    sntracker.draw_content_label('{:.2f}')

    for i_cell in range(len(stds)):
        if stds[i_cell] is None:
            continue
        sntracker.setcontent(i_cell, stds[i_cell])

    sntracker.draw()
    sntracker.save(output_directory)

    sntracker = sn.tracker(name + "_chi2", True)
    sntracker.draw_cellid_label()
    sntracker.draw_content_label('{:.2f}')

    for i_cell in range(len(stds)):
        if chis[i_cell] is None:
            continue
        sntracker.setcontent(i_cell, chis[i_cell])

    sntracker.draw()
    sntracker.save(output_directory)

    return means


def plot_drift_time_vs_voltage(directory: str, output_directory):
    files = ["red_608_output.root", "red_617_output.root", "red_618_output.root", "red_619_output.root",
             "red_623_output.root", "red_635_output.root", "red_638_output.root", "red_640_output.root",
             "red_641_output.root", "red_642_output.root", "red_645_output.root", "red_652_output.root",
             "red_653_output.root", "red_654_output.root", "red_660_output.root", "red_668_output.root",
             "red_669_output.root", "red_670_output.root", "red_671_output.root"]
    times = [60, 10, 10, 10,
             15, 3.66667, 30, 5,
             5, 5, 15, 10,
             10, 10, 15, 10,
             10, 10, 10]
    HVs = [1600, 1600, 1600, 1600,
           1600, 1650, 1612, 1650,
           1675, 1700, 1625, 1650,
           1650, 1650, 1675, 1600,
           1625, 1650, 1675]
    Areas = [0, 3, 3, 3,
             4, 5, 6, 6,
             6, 6, 7, 1,
             1, 1, 2, 4,
             4, 4, 4]
    to_do = [1, 1, 0, 0,
             0, 1, 1, 1,
             1, 1, 1, 0,
             0, 1, 1, 1,
             1, 1, 1]
    is_first = [True for i in range(8)]
    temps = [-2.5, -2.5, -2.5, -5, 2.5, 5, 0, 0]
    markers = ['s', 'o']
    data = [[], []]

    for i_file, file in enumerate(files):
        if to_do[i_file] == 0:
            continue
        drifts = [[], []]
        root_file = ROOT.TFile(directory + "/" + file, "READ")
        event_tree = root_file.event_tree
        for event in event_tree:
            if len(event.tracker_time_anode) > 50:
                continue
            '''if len(event.calo_time) == 0:
                continue
            ht = False
            for j in range(len(event.calo_time)):
                if bool(event.calo_high_t[j]):
                    ht = True
            if not ht:
                continue'''

            for i in range(len(event.tracker_time_anode)):
                if event.tracker_time_top_cathode[i] != is_none and event.tracker_time_bottom_cathode[i] != is_none and event.tracker_time_anode[i] != is_none:
                    time = event.tracker_time_anode[i]
                    drift = ((event.tracker_time_top_cathode[i] - time) + (event.tracker_time_bottom_cathode[i] - time))*1e6

                    side, row, layer = cell_id(event.tracker_cell_num[i])
                    if layer in [0, 8] or row in [0, 13, 14, 27, 28, 41, 42, 55, 56, 68, 71, 84, 85, 98, 99, 112]:
                        drifts[0].append(drift)
                    else:
                        drifts[1].append(drift)

        lower, upper = 0, 100
        n_bins = 50
        names = ['Edges', 'Centres']
        fig = plt.figure(figsize=figsize, facecolor='white')
        for i in range(len(drifts)):
            if len(drifts[i]) == 0:
                continue

            mean = np.average(drifts[i])
            std = np.std(drifts[i])

            freq, bin_edges = np.histogram(drifts[i], n_bins, range=(lower, upper))
            width = bin_edges[-1] - bin_edges[-2]
            bin_centres = bin_edges[:-1] + width / 2

            fit = ROOT.TF1("fit", "[0]*TMath::Gaus(x, [1], [2]) + [3]", lower, upper)
            hist = ROOT.TH1D("", "", n_bins, lower, upper)
            for j in range(len(drifts[i])):
                hist.Fill(drifts[i][j])
            fit.SetParLimits(0, 0, (np.max(freq) + np.max(freq)*0.1))
            if bin_centres[np.argmax(freq)] - 10 < 0 or bin_centres[np.argmax(freq)] + 10 > 100:
                fit.SetParLimits(1, 0, 100)
            else:
                fit.SetParLimits(1, bin_centres[np.argmax(freq)] - 10, bin_centres[np.argmax(freq)] + 10)
            fit.SetParLimits(2, 0, 50)
            fit.SetParameters(np.max(freq), bin_centres[np.argmax(freq)], 1, 1)
            hist.Fit("fit", "0Q")
            A = fit.GetParameter(0)
            mu = fit.GetParameter(1)
            sig = fit.GetParameter(2)
            B = fit.GetParameter(3)
            x = np.linspace(lower, upper, 500)

            if fit.GetNDF() == 0:
                chi = None
            else:
                chi = fit.GetChisquare() / fit.GetNDF()

            plt.bar(bin_centres, freq, width=width, alpha=0.5, label=names[i])
            plt.plot(x, gaussian_B(x, A, mu, sig, B), 'C{}-'.format(i))
            plt.axvline(mu, ls='--', color='C{}'.format(i))

            data[i].append((HVs[i_file], mu, sig, Areas[i_file]))

            del hist
            del fit

        plt.title('Area {} run {}'.format(Areas[i_file], file.split("_")[1]))
        plt.xlabel("Plasma Propagation Time /µs")
        plt.ylabel("Counts")
        plt.legend(loc='best')
        plt.xlim(0, 100)
        plt.tight_layout()
        plt.savefig(output_directory + "/run_{}_all_cell_drift_times.pdf".format(file.split("_")[1]))
        plt.close()

        '''if is_first[Areas[i_file]]:
                    plt.plot(HVs[i_file], mu, color='C{}'.format(Areas[i_file]), marker='.', label='Area {}'.format(Areas[i_file]))
                    is_first[Areas[i_file]] = False
                else:
                    plt.plot(HVs[i_file], mu, color='C{}'.format(Areas[i_file]), marker='.')

                del hist
                del fit'''

    fig = plt.figure(figsize=figsize, facecolor='white')
    for i in range(len(data)):
        for j in range(len(data[i])):
            if is_first[data[i][j][3]]:
                is_first[data[i][j][3]] = False
                plt.errorbar(data[i][j][0] + temps[data[i][j][3]], data[i][j][1], yerr=data[i][j][2],
                             fmt='C{}{}'.format(data[i][j][3], markers[i]), label='Area {}'.format(data[i][j][3]),
                             markersize=3)
            else:
                plt.errorbar(data[i][j][0] + temps[data[i][j][3]], data[i][j][1], yerr=data[i][j][2],
                             fmt='C{}{}'.format(data[i][j][3], markers[i]), markersize=3)
    plt.plot(0, 0, 'ks', label='Edges')
    plt.plot(0, 0, 'ko', label='Centres')
    plt.xlabel("Voltage /V")
    plt.ylabel("Plasma Propagation Time /µs")
    # plt.title("Run {} All Cells".format(run_num))
    plt.legend(bbox_to_anchor=(1.05, 1))
    plt.xlim(1575, 1725)
    plt.ylim(0, 100)
    plt.tight_layout()
    plt.savefig(output_directory + "/all_drift_times_vs_voltage.pdf")
    plt.close()


def plot_ppt_vs_dt(directory, output_directory):
    files = ["red_608_output.root", "red_617_output.root", "red_618_output.root", "red_619_output.root",
             "red_623_output.root", "red_635_output.root", "red_638_output.root", "red_640_output.root",
             "red_641_output.root", "red_642_output.root", "red_645_output.root", "red_652_output.root",
             "red_653_output.root", "red_654_output.root", "red_660_output.root", "red_668_output.root",
             "red_669_output.root", "red_670_output.root", "red_671_output.root", "red_612_output.root"]
    times = [60, 10, 10, 10,
             15, 3.66667, 30, 5,
             5, 5, 15, 10,
             10, 10, 15, 10,
             10, 10, 10, 60]
    HVs = [1600, 1600, 1600, 1600,
           1600, 1650, 1612, 1650,
           1675, 1700, 1625, 1650,
           1650, 1650, 1675, 1600,
           1625, 1650, 1675, 1600]
    Areas = [0, 3, 3, 3,
             4, 5, 6, 6,
             6, 6, 7, 1,
             1, 1, 2, 4,
             4, 4, 4, 3]
    to_do = [1, 0, 0, 0,
             0, 0, 1, 0,
             0, 0, 0, 0,
             0, 0, 0, 0,
             0, 0, 0, 1]

    is_first = [True for i in range(8)]
    temps = [-2.5, -2.5, -2.5, -5, 2.5, 5, 0, 0]
    markers = ['s', 'o']

    fig = plt.figure(figsize=(6,3), facecolor='white')
    for i_file, file in enumerate(files):
        dts = [0 for i in range(2034)]
        drifts = [[[], []], [[], []]]
        if to_do[i_file] == 0:
            continue
        root_file = ROOT.TFile(directory + "/" + file, "READ")
        event_tree = root_file.event_tree
        for event in event_tree:
            if len(event.tracker_time_anode) > 50:
                continue
            '''if len(event.calo_time) != 0:
                continue
                ht = False
                for j in range(len(event.calo_time)):
                    if bool(event.calo_high_t[j]):
                        ht = True
                if not ht:
                    continue'''

            for i in range(len(event.tracker_time_anode)):
                if event.tracker_time_top_cathode[i] != is_none and event.tracker_time_bottom_cathode[i] != is_none and event.tracker_time_anode[i] != is_none:
                    time = event.tracker_time_anode[i]
                    dt = time*1e6 - dts[event.tracker_cell_num[i]]
                    if dts[event.tracker_cell_num[i]] == 0:
                        dts[event.tracker_cell_num[i]] = event.tracker_time_anode[i]*1e6
                        continue
                    dts[event.tracker_cell_num[i]] = time * 1e6
                    drift = ((event.tracker_time_top_cathode[i] - time) + (event.tracker_time_bottom_cathode[i] - time)) * 1e6
                    if dt > 10000000:
                        continue
                    side, row, layer = cell_id(event.tracker_cell_num[i])
                    if layer in [0, 8] or row in [0, 13, 14, 27, 28, 41, 42, 55, 56, 68, 71, 84, 85, 98, 99, 112]:
                        drifts[0][0].append(dt)
                        drifts[0][1].append(drift)
                    else:
                        drifts[1][0].append(dt)
                        drifts[1][1].append(drift)

        if is_first[Areas[i_file]]:
            plt.plot(np.array(drifts[0][0]) / 1000, drifts[0][1], 'C{}s'.format(Areas[i_file]),
                     label='Area {}'.format(Areas[i_file]), markersize=2)
            plt.plot(np.array(drifts[1][0]) / 1000, drifts[1][1], 'C{}o'.format(Areas[i_file]), markersize=2)
        else:
            plt.plot(np.array(drifts[0][0]) / 1000, drifts[0][1], 'C{}s'.format(Areas[i_file]), markersize=2)
            plt.plot(np.array(drifts[1][0]) / 1000, drifts[1][1], 'C{}o'.format(Areas[i_file]), markersize=2)

    plt.plot(-1, -1, 'ks', label='Edges')
    plt.plot(-1, -1, 'ko', label='Centres')
    plt.legend(loc='upper right')
    plt.xlim(0, 50)
    plt.ylim(40, 100)
    plt.xlabel(r"$\Delta t$ /ms")
    plt.ylabel("Plasma Propagation Time /µs")
    plt.tight_layout()
    plt.savefig(output_directory + "/ppt_vs_dt.pdf")


def plot_r56_vs_r12(directory, file, output_directory, run_time):
    root_file = ROOT.TFile(directory + "/" + file, "READ")
    run_num = file.split("_")[1]
    event_tree = root_file.event_tree
    drifts = [[], [], []]
    drifts_56 = []
    drifts_12 = []
    drifts_34 = []
    first = [[], [], []]
    second = [[], [], []]
    for event in event_tree:
        if len(event.tracker_time_anode) > 50:
            continue

        for i in range(len(event.tracker_time_anode)):
            if event.tracker_time_top_cathode[i] != is_none and event.tracker_time_bottom_cathode[i] != is_none and \
                    event.tracker_time_anode[i] != is_none:
                time = event.tracker_time_anode[i]

                t1 = (event.tracker_time_anode_first_lt[i] - time) * 1e6
                t2 = (event.tracker_time_anode_first_ht[i] - time) * 1e6
                t3 = (event.tracker_time_anode_second_lt[i] - time) * 1e6
                t4 = (event.tracker_time_anode_second_ht[i] - time) * 1e6
                t5 = (event.tracker_time_bottom_cathode[i] - time) * 1e6
                t6 = (event.tracker_time_top_cathode[i] - time) * 1e6

                drift_56 = t5 + t6
                drifts_56.append(drift_56)

                if t6 < t5:
                    first[0].append(t6)
                    second[0].append(t5)
                else:
                    second[0].append(t6)
                    first[0].append(t5)

                drifts[0].append(drift_56)
                if not (0 in [t1, t2]):
                    drift_12 = t1 + t2
                    drifts_12.append(drift_12)
                    drifts[1].append(drift_12)
                else:
                    drifts_12.append(0)
                    drifts[1].append(0)
                if not (0 in [t3, t4]):
                    drift_34 = t3 + t4
                    drifts_34.append(drift_34)
                    drifts[2].append(drift_34)
                else:
                    drifts_34.append(0)
                    drifts[2].append(0)

                first[1].append(t1)
                first[2].append(t3)
                second[1].append(t2)
                second[2].append(t4)

    n_bins = 100
    lower, upper = 30, 80
    freq, bin_edges = np.histogram(drifts_56, n_bins, range=(lower, upper))
    width = bin_edges[-1] - bin_edges[-2]
    bin_centres = bin_edges[:-1] + width / 2
    fig = plt.figure(figsize=figsize, facecolor='white')
    plt.bar(bin_centres, freq / run_time, width=width, color='C0')
    plt.xlabel("Plasma Propagation Time /µs")
    plt.ylabel("Counts per Minute")
    plt.title("Run: " + run_num + " PPT - 56")
    plt.tight_layout()
    plt.savefig(output_directory + "/plots/run_" + run_num + "/" + run_num + '_all_drift_56.pdf')
    plt.close()

    freq, bin_edges = np.histogram(drifts_12, n_bins, range=(lower, upper))
    width = bin_edges[-1] - bin_edges[-2]
    bin_centres = bin_edges[:-1] + width / 2

    fig = plt.figure(figsize=figsize, facecolor='white')
    plt.bar(bin_centres, freq / run_time, width=width, color='C1')
    plt.title("Run: " + run_num + " PPT - 12")
    plt.xlabel("Plasma Propagation Time /µs")
    plt.ylabel("Counts per Minute")
    plt.xlim(lower, upper)
    plt.tight_layout()
    plt.savefig(output_directory + "/plots/run_" + run_num + "/" + run_num + '_all_drift_12.pdf')
    plt.close()

    freq, bin_edges = np.histogram(drifts_34, n_bins, range=(lower, upper))
    width = bin_edges[-1] - bin_edges[-2]
    bin_centres = bin_edges[:-1] + width / 2

    fig = plt.figure(figsize=figsize, facecolor='white')
    plt.bar(bin_centres, freq / run_time, width=width, color='C2')
    plt.title("Run: " + run_num + " PPT - 34")
    plt.xlabel("Plasma Propagation Time /µs")
    plt.ylabel("Counts per Minute")
    plt.xlim(lower, upper)
    plt.tight_layout()
    plt.savefig(output_directory + "/plots/run_" + run_num + "/" + run_num + '_all_drift_34.pdf')
    plt.close()

    fig, ax = plt.subplots(figsize=figsize, facecolor='white')
    counts, xedges, yedges, im = ax.hist2d(drifts[0], drifts[1], bins=(20, 20), range=[[40, 80], [40, 80]], cmin=1,
                                           cmap=plt.get_cmap('Reds'))
    plt.title("Run: " + run_num + " PPT - 56 vs 12")
    plt.xlabel("t56 PPT /µs")
    plt.ylabel("t12 PPT /µs")
    plt.xlim(40, 80)
    plt.ylim(40, 80)
    fig.colorbar(im, ax=ax)
    plt.tight_layout()
    plt.savefig(output_directory + "/plots/run_" + run_num + "/" + run_num + '_all_drift_56_12.pdf')
    plt.close()

    fig, ax = plt.subplots(figsize=figsize, facecolor='white')
    counts, xedges, yedges, im = ax.hist2d(drifts[0], drifts[2], bins=(20, 20), range=[[40, 80], [40, 80]], cmin=1,
                                           cmap=plt.get_cmap('Blues'))
    plt.title("Run: " + run_num + " PPT - 56 vs 34")
    plt.xlabel("t56 PPT /µs")
    plt.ylabel("t34 PPT /µs")
    plt.xlim(40, 80)
    plt.ylim(40, 80)
    fig.colorbar(im, ax=ax)
    plt.tight_layout()
    plt.savefig(output_directory + "/plots/run_" + run_num + "/" + run_num + '_all_drift_56_34.pdf')
    plt.close()

    fig, ax = plt.subplots(figsize=figsize, facecolor='white')
    counts, xedges, yedges, im = ax.hist2d(first[0], first[1], bins=(20, 20), range=[[0, 100], [0, 100]], cmin=1,
                                           cmap=plt.get_cmap('Blues'))
    plt.title("Run: " + run_num + " 1st Cathode vs LT Anode")
    plt.xlabel("First Cathode time /µs")
    plt.ylabel("t1 /µs")
    plt.xlim(0, 100)
    plt.ylim(0, 100)
    fig.colorbar(im, ax=ax)
    plt.tight_layout()
    plt.savefig(output_directory + "/plots/run_" + run_num + "/" + run_num + '_first_t1.pdf')
    plt.close()

    fig, ax = plt.subplots(figsize=figsize, facecolor='white')
    counts, xedges, yedges, im = ax.hist2d(first[0], first[2], bins=(20, 20), range=[[0, 100], [0, 100]], cmin=1,
                                           cmap=plt.get_cmap('Reds'))
    plt.title("Run: " + run_num + " 1st Cathode vs HT Anode")
    plt.xlabel("First Cathode time /µs")
    plt.ylabel("t3 /µs")
    plt.xlim(0, 100)
    plt.ylim(0, 100)
    fig.colorbar(im, ax=ax)
    plt.tight_layout()
    plt.savefig(output_directory + "/plots/run_" + run_num + "/" + run_num + '_first_t3.pdf')
    plt.close()

    fig, ax = plt.subplots(figsize=figsize, facecolor='white')
    counts, xedges, yedges, im = ax.hist2d(second[0], second[1], bins=(20, 20), range=[[0, 100], [0, 100]], cmin=1,
                                           cmap=plt.get_cmap('Greens'))
    plt.title("Run: " + run_num + " 2nd Cathode vs LT Anode")
    plt.xlabel("Second Cathode time /µs")
    plt.ylabel("t2 /µs")
    plt.xlim(0, 100)
    plt.ylim(0, 100)
    fig.colorbar(im, ax=ax)
    plt.tight_layout()
    plt.savefig(output_directory + "/plots/run_" + run_num + "/" + run_num + '_second_t2.pdf')
    plt.close()

    fig, ax = plt.subplots(figsize=figsize, facecolor='white')
    counts, xedges, yedges, im = ax.hist2d(second[0], second[2], bins=(20, 20), range=[[0, 100], [0, 100]], cmin=1,
                                           cmap=plt.get_cmap('Purples'))
    plt.title("Run: " + run_num + " 2nd Cathode vs HT Anode")
    plt.xlabel("Second Cathode time /µs")
    plt.ylabel("t4 /µs")
    plt.xlim(0, 100)
    plt.ylim(0, 100)
    fig.colorbar(im, ax=ax)
    plt.tight_layout()
    plt.savefig(output_directory + "/plots/run_" + run_num + "/" + run_num + '_second_t4.pdf')
    plt.close()

####################################
# New functions
####################################


def parse_root_file(file_name):
    root_file = ROOT.TFile(file_name, "READ")
    tree = root_file.RED

    # Build a dictionary full of the data you want to keep
    data = {}

    n_events = tree.GetEntries()
    counter = 0
    print(f">>> Events in tree: {n_events}")
    print(">>> Parsing ...")
    print(">")
    for event in tree:
        ################################################################################################################
        #
        #  Access the branches of the root file
        #  To know what is in the TTree either use:
        #               tree.Print()
        #               OR
        #               branches = tree.GetListOfBranches()
        #               print([branch.GetName() for branch in branches])
        #
        ################################################################################################################
        run_id = event.run_id
        event_id = event.event_id
        # event_clock = event.event_clock                                       # int The clock of the events seconds
        # event_ticks = event.event_ticks                                       # Event times - TODO: not yet implemented
        n_calo_hits = int(event.nb_calo_hits)                                   # int Number of calorimeter hits
        # calo_hit_id = list(event.calo_hit_id)                                 # list(int) List of the calo hit IDs
        calo_om_side_id = list(event.calo_om_side_id)                           # list(int) List of OM sides (1: French, 0: Italian)
        calo_om_wall_id = list(event.calo_om_wall_id)                           # list(int) List of OM walls (0,1) see url:https://nemo.lpc-caen.in2p3.fr/wiki/NEMO/SuperNEMO/Calorimeter
        calo_om_column_id = list(event.calo_om_column_id)                       # list(int) List of OM columns
        calo_om_row_id = list(event.calo_om_row_id)                             # list(int) List of OM rows
        # calo_om_ref_id = list(event.calo_om_ref_id)                           # Don't know?
        # calo_clock = list(event.calo_clock)                                   # Don't know?
        calo_ticks = list(event.calo_ticks)                                     # list(int) Calo TDC (6.25ns)
        calo_lto = list(event.calo_lto)                                         # list(bool) Calo low threshold flag
        calo_ht = list(event.calo_ht)                                           # list(bool) Calo high threshold flag
        # calo_fcr = list(event.calo_fcr)                                       # list(int) Dont' know
        # calo_lt_trigger_counter = list(event.calo_lt_trigger_counter)         # Don't know
        # calo_lt_time_counter = list(event.calo_lt_time_counter)               # Don't know
        calo_fwmeas_baseline = list(event.calo_fwmeas_baseline)                 # list(int) Calo baseline LSB: ADC unit/16
        calo_fwmeas_peak_amplitude = list(event.calo_fwmeas_peak_amplitude)     # list(int) Calo amplitude LSB: ADC unit/8
        calo_fwmeas_peak_cell = list(event.calo_fwmeas_peak_cell)               # list(int) Calo peak cell TDC: 0-1023
        calo_fwmeas_charge = list(event.calo_fwmeas_charge)                     # list(int) Calo charge NO IDEA OF UNITS - do not trust this
        calo_fwmeas_rising_cell = list(event.calo_fwmeas_rising_cell)           # list(int) Calo rising cell (right edge) TDC unit/256
        calo_fwmeas_falling_cell = list(event.calo_fwmeas_falling_cell)         # list(int) Calo falling cell (left edge) TDC unit/256

        calo_waveform = [list(event.calo_waveform)[1024 * i:1024 * (i + 1)] for i in range(n_calo_hits)]  # list(list(int)) Calo waveform

        nb_tracker_hits = event.nb_tracker_hits                                 # int Number of tracker cells that have at least one hit
        # tracker_hit_id = list(event.tracker_hit_id)                           # list(int) TR hit ids
        tracker_cell_side_id = list(event.tracker_cell_side_id)                 # list(int) TR side (1: French, 0: Italian)
        tracker_cell_row_id = list(event.tracker_cell_row_id)                   # list(int) TR row (0-112)
        tracker_cell_layer_id = list(event.tracker_cell_layer_id)               # list(int) TR layer (0-8)
        # tracker_clock = list(event.tracker_clock)                             # list(int) Don't know

        tracker_anode_R0_ticks = [list(event.tracker_anode_R0_ticks)[i * MAX_GG_TIMES:MAX_GG_TIMES * (i + 1)] for i in
                                  range(nb_tracker_hits)]  # list(int) Anode time LSB: 12.5 ns
        tracker_anode_R1_ticks = [list(event.tracker_anode_R1_ticks)[i * MAX_GG_TIMES:MAX_GG_TIMES * (i + 1)] for i in
                                  range(nb_tracker_hits)]  # list(int) Anode: 1st low threshold LSB: 12.5 ns
        tracker_anode_R2_ticks = [list(event.tracker_anode_R2_ticks)[i * MAX_GG_TIMES:MAX_GG_TIMES * (i + 1)] for i in
                                  range(nb_tracker_hits)]  # list(int) Anode: 2nd low threshold LSB: 12.5 ns
        tracker_anode_R3_ticks = [list(event.tracker_anode_R3_ticks)[i * MAX_GG_TIMES:MAX_GG_TIMES * (i + 1)] for i in
                                  range(nb_tracker_hits)]  # list(int) Anode: 1st high threshold LSB: 12.5 ns
        tracker_anode_R4_ticks = [list(event.tracker_anode_R4_ticks)[i * MAX_GG_TIMES:MAX_GG_TIMES * (i + 1)] for i in
                                  range(nb_tracker_hits)]  # list(int) Anode: 2st high threshold LSB: 12.5 ns
        tracker_bottom_cathode_R5_ticks = [list(event.tracker_bottom_cathode_R5_ticks)[i * MAX_GG_TIMES:MAX_GG_TIMES * (i + 1)] for i in
                                           range(nb_tracker_hits)]  # list(int) Cathode bottom LSB: 12.5 ns
        tracker_top_cathode_R6_ticks = [list(event.tracker_top_cathode_R6_ticks)[i * MAX_GG_TIMES:MAX_GG_TIMES * (i + 1)] for i in
                                        range(nb_tracker_hits)]  # list(int) Cathode top LSB: 12.5 ns

        ################################################################################################################
        #  Process data here either store in a container like in the data dictionary
        #  OR process directly here
        ################################################################################################################

        ################
        # Quality Cuts #
        ################
        if n_calo_hits == 0 or nb_tracker_hits == 0 or nb_tracker_hits > 50:
            # Cut out noisy data
            continue

        calo_om_nums = [None for i in range(n_calo_hits)]
        for i_om in range(n_calo_hits):
            if calo_om_row_id[i_om] == -1:
                # GVETO OM
                calo_om_num = 520 + 128 + calo_om_side_id[i_om] * 32 + calo_om_wall_id[i_om] * 16 + calo_om_column_id[i_om]
            elif calo_om_wall_id[i_om] == -1:
                # MAIN WALL OM
                calo_om_num = calo_om_side_id[i_om] * 20 * 13 + calo_om_column_id[i_om] * 13 + calo_om_row_id[i_om]
            else:
                # XWALL OM
                calo_om_num = 520 + calo_om_side_id[i_om] * 64 + calo_om_wall_id[i_om] * 32 + calo_om_column_id[i_om] * 16 + calo_om_row_id[i_om]
            calo_om_nums[i_om] = calo_om_num

        try:
            # Take the first calorimeter time
            calo_time = min(calo_ticks) * c_tdc2sec
        except ValueError:
            # If there are no calo hits (its likely noise) set the calo time to 0
            calo_time = 0

        tracker_cell_nums = [None for i in range(nb_tracker_hits)]
        tracker_ppts = [None for i in range(nb_tracker_hits)]
        tracker_cell_times = [[] for i in range(nb_tracker_hits)]
        tracker_cell_timestamps = [[] for i in range(nb_tracker_hits)]

        for i_cell in range(nb_tracker_hits):
            tracker_cell_num = tracker_cell_side_id[i_cell] * 113 * 9 + tracker_cell_row_id[i_cell] * 9 + tracker_cell_layer_id[i_cell]
            # Note: There are 10 entries for each of the cell R0-6. I have chosen to use the 0th entry
            # >>> If the value is -9223372036854775808 this is the default value
            # >>> It is when there is a value for one RX but not another RY for a specific cell - ignore these values

            # TDC RX
            rs = [tracker_anode_R0_ticks[i_cell][0],
                  tracker_anode_R1_ticks[i_cell][0], tracker_anode_R2_ticks[i_cell][0],
                  tracker_anode_R3_ticks[i_cell][0], tracker_anode_R4_ticks[i_cell][0],
                  tracker_bottom_cathode_R5_ticks[i_cell][0], tracker_top_cathode_R6_ticks[i_cell][0]]
            tracker_cell_timestamps[i_cell] = rs
            # Time (RX - calo hit time) in µs
            ts = []
            for r in rs:
                if r == -9223372036854775808:
                    ts.append(None)
                elif r == 0:
                    ts.append(None)
                else:
                    ts.append((r*t_tdc2sec - calo_time)*1E6)
            tracker_cell_times[i_cell] = ts
            tracker_cell_nums[i_cell] = tracker_cell_num
            try:
                tracker_ppts[i_cell] = ts[5] + ts[6]
            except TypeError:
                pass

        # Fill the dictionary
        data[event_id] = {
            "cells": {cell: {
                "rs": tracker_cell_timestamps[index],
                "ts": tracker_cell_times[index],
                "ppts": tracker_ppts[index]
            } for index, cell in enumerate(tracker_cell_nums)},
            "oms": {om: {
                "amplitudes": -1 * calo_fwmeas_peak_amplitude[index] * adc2mv/8
            } for index, om in enumerate(calo_om_nums)}
        }

        '''if event_id == 0:
            for index in range(len(tracker_cell_timestamps)):
                print(tracker_cell_nums[index], tracker_cell_timestamps[index])'''

        # Progress
        new_counter = (event_id + 1) / n_events * 100 // 10
        if new_counter != counter:
            counter = new_counter
            print("> Progress: {:.1f}".format(event_id / n_events * 100) + f"%. Events {len(data.keys())}")
    return data


def plot_3d_event(filename, ppts):
    root_file = ROOT.TFile(filename, "READ")
    event_tree = root_file.event_tree
    D = 2
    n_event = 0
    plot_event = False
    for event in event_tree:
        om_hit = [False for i in range(712)]
        anode_times = [None for i in range(2034)]
        top_cathode_times = [None for i in range(2034)]
        bot_cathode_times = [None for i in range(2034)]
        pos = [None for i in range(2034)]
        if len(event.tracker_time_anode) > 50:
            continue
        calo_time = None
        ht = False
        if len(event.calo_time) != 0:
            for j in range(len(event.calo_time)):
                if bool(event.calo_high_t[j]):
                    ht = True
                    temp_calo_time = event.calo_time[j]
                    if calo_time is not None and temp_calo_time < calo_time:
                        calo_time = temp_calo_time
                    elif calo_time is None:
                        calo_time = temp_calo_time
                    om_hit[event.calo_om_num[j]] = True
                    # plot_event = True

        for i in range(len(event.tracker_time_anode)):
            if event.tracker_time_top_cathode[i] != is_none and event.tracker_time_bottom_cathode[i] != is_none and \
                    event.tracker_time_anode[i] != is_none:
                time = event.tracker_time_anode[i]
                if calo_time is not None:
                    anode_times.append((event.tracker_time_anode[i] - calo_time) * 1e6)
                t5 = (event.tracker_time_bottom_cathode[i] - time) * 1e6
                t6 = (event.tracker_time_top_cathode[i] - time) * 1e6
                top_cathode_times[event.tracker_cell_num[i]] = t6
                bot_cathode_times[event.tracker_cell_num[i]] = t5

                if ppts[event.tracker_cell_num[i]] is None:
                    continue
                if t5 < t6:
                    d = D * t5 / ppts[event.tracker_cell_num[i]]
                else:
                    d = D*(1 - t6 / ppts[event.tracker_cell_num[i]])
                pos[event.tracker_cell_num[i]] = d

        if plot_event:
            sndemo = sn.demonstrator("event_" + str(n_event))
            n_calo = 0
            n_tracker = 0
            for om in range(712):
                if om_hit[om]:
                    n_calo += 1
                    sndemo.setomcontent(om, 1)
            for cell in range(2034):
                if pos[cell] is not None:
                    n_tracker += 1
                    sndemo.setggcontent(cell, 1)
            if n_calo == 2 and n_tracker > 5:
                try:
                    sndemo.draw_top()
                    sndemo.save("/Users/williamquinn/Desktop/read_red")
                except ValueError:
                    print("sndisplay problem")
            del sndemo

        if n_event == 88:
            with open("/Users/williamquinn/Desktop/read_red/event_88.csv", "w") as outfile:
                for om in range(712):
                    if om_hit[om]:
                        outfile.write("OM," + str(om) + "\n")
                for cell in range(2034):
                    if pos[cell] is not None:
                        outfile.write("GG," + str(cell) + "," + str(pos[cell]) + "\n")

        n_event += 1
    root_file.Close()


def plot_all_times(file_name, run_time):
    run = file_name.split("/")[-1].split("_")[1]
    file = ROOT.TFile(file_name, "READ")
    tree = file.event_tree
    r0s, r1s, r2s, r3s, r4s, r5s, r6s = [], [], [], [], [], [], []

    for event in tree:
        if len(event.tracker_time_anode) > 50:
            continue
        if len(list(event.calo_time)) == 0:
            continue
        for i in range(len(event.tracker_time_anode)):
            if event.tracker_time_top_cathode[i] != is_none and event.tracker_time_bottom_cathode[i] != is_none and \
                    event.tracker_time_anode[i] != is_none:
                time = np.amin(list(event.calo_time)) * 1e6
                r0 = event.tracker_time_anode[i] * 1e6 - time
                r1 = event.tracker_time_anode_first_lt[i] * 1e6 - time
                r2 = event.tracker_time_anode_second_lt[i] * 1e6 - time
                r3 = event.tracker_time_anode_first_ht[i] * 1e6 - time
                r4 = event.tracker_time_anode_second_ht[i] * 1e6 - time
                r5 = event.tracker_time_bottom_cathode[i] * 1e6 - time
                r6 = event.tracker_time_top_cathode[i] * 1e6 - time

                r0s.append(r0)
                r1s.append(r1)
                r2s.append(r2)
                r3s.append(r3)
                r4s.append(r4)
                r5s.append(r5)
                r6s.append(r6)
    rs = [r0s, r1s, r2s, r3s, r4s, r5s, r6s]
    plt.figure(figsize=figsize)
    lower, higher = -20, 100
    lw = 1
    for index, r in enumerate(rs):
        freq, bin_edges = np.histogram(r, range=(lower, higher), bins=120)
        width = bin_edges[-1] - bin_edges[-2]
        bin_centres = bin_edges[:-1] + width/2
        plt.plot(bin_centres, freq/run_time, "C{}".format(index), label='R{}'.format(index), linewidth=lw)

    plt.xlabel(r"Time /$\mu$s")
    plt.xlim(-40, higher)
    plt.title("All Cells Timestamps R0-6")
    plt.legend(loc='best')
    plt.ylabel("Counts per minute")
    plt.yscale("log")
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/read_red/{}_all_times.pdf".format(run))

    file.Close()


def plot_anode_time_comps(file_name, run_time):
    run = file_name.split("/")[-1].split("_")[1]
    file = ROOT.TFile(file_name, "READ")
    tree = file.event_tree
    r0s, r1s, r2s, r3s, r4s, r5s, r6s = [], [], [], [], [], [], []

    for event in tree:
        if len(event.tracker_time_anode) > 50:
            continue
        if len(list(event.calo_time)) == 0:
            continue
        for i in range(len(event.tracker_time_anode)):
            if event.tracker_time_top_cathode[i] != is_none and event.tracker_time_bottom_cathode[i] != is_none and \
                    event.tracker_time_anode[i] != is_none:
                time = np.amin(list(event.calo_time)) * 1e6
                r0 = event.tracker_time_anode[i] * 1e6 - time
                r1 = event.tracker_time_anode_first_lt[i] * 1e6 - time
                r2 = event.tracker_time_anode_second_lt[i] * 1e6 - time
                r3 = event.tracker_time_anode_first_ht[i] * 1e6 - time
                r4 = event.tracker_time_anode_second_ht[i] * 1e6 - time
                r5 = event.tracker_time_bottom_cathode[i] * 1e6 - time
                r6 = event.tracker_time_top_cathode[i] * 1e6 - time

                r0s.append(r0)
                r1s.append(r1)
                r2s.append(r2)
                r3s.append(r3)
                r4s.append(r4)
                r5s.append(r5)
                r6s.append(r6)
    file.Close()

    lower, higher = 0, 100
    n = 50

    # t1 vs t2
    plt.figure(figsize=(3, 2))
    plt.hist2d(r1s, r2s, bins=(n, n), range=[[lower, higher], [lower, higher]], cmap='Blues')
    plt.xlabel(r"Timestamp t1 /$\mu$s")
    plt.ylabel(r"Timestamp t2 /$\mu$s")
    plt.colorbar()
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/read_red/{}_t1_v_t2_times.pdf".format(run))

    # t1 vs t3
    plt.figure(figsize=(3, 2))
    plt.hist2d(r1s, r3s, bins=(n, n), range=[[lower, higher], [lower, higher]], cmap='Greens')
    plt.xlabel(r"Timestamp t1 /$\mu$s")
    plt.ylabel(r"Timestamp t3 /$\mu$s")
    plt.colorbar()
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/read_red/{}_t1_v_t3_times.pdf".format(run))

    # t2 vs t4
    plt.figure(figsize=(3, 2))
    plt.hist2d(r2s, r4s, bins=(n, n), range=[[lower, higher], [lower, higher]], cmap='Purples')
    plt.xlabel(r"Timestamp t2 /$\mu$s")
    plt.ylabel(r"Timestamp t4 /$\mu$s")
    plt.colorbar()
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/read_red/{}_t2_v_t4_times.pdf".format(run))

    # t3 vs t4
    plt.figure(figsize=(3, 2))
    plt.hist2d(r3s, r4s, bins=(n, n), range=[[lower, higher], [lower, higher]], cmap='Reds')
    plt.xlabel(r"Timestamp t3 /$\mu$s")
    plt.ylabel(r"Timestamp t4 /$\mu$s")
    plt.colorbar()
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/read_red/{}_t3_v_t4_times.pdf".format(run))

    plt.figure(figsize=figsize)
    plt.hist2d(r5s, r6s, bins=(n, n), range=[[lower, higher], [lower, higher]], cmap='Oranges')
    plt.xlabel(r"Timestamp t5 /$\mu$s")
    plt.ylabel(r"Timestamp t6 /$\mu$s")
    plt.colorbar()
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/read_red/{}_t5_v_t6_times.pdf".format(run))


def plot_ppt_all(file_name, run_time):
    run = file_name.split("/")[-1].split("_")[1]
    file = ROOT.TFile(file_name, "READ")
    tree = file.event_tree
    ppts = []

    for event in tree:
        if len(event.tracker_time_anode) > 50:
            continue
        if len(list(event.calo_time)) == 0:
            continue
        for i in range(len(event.tracker_time_anode)):
            if event.tracker_time_top_cathode[i] != is_none and event.tracker_time_bottom_cathode[i] != is_none and \
                    event.tracker_time_anode[i] != is_none:
                time = np.amin(list(event.calo_time)) * 1e6
                r0 = event.tracker_time_anode[i] * 1e6 - time
                r1 = event.tracker_time_anode_first_lt[i] * 1e6 - time
                r2 = event.tracker_time_anode_second_lt[i] * 1e6 - time
                r3 = event.tracker_time_anode_first_ht[i] * 1e6 - time
                r4 = event.tracker_time_anode_second_ht[i] * 1e6 - time
                r5 = event.tracker_time_bottom_cathode[i] * 1e6 - time
                r6 = event.tracker_time_top_cathode[i] * 1e6 - time
                ppts.append(r5 + r6)
    file.Close()

    plt.figure(figsize=figsize)
    freq, bin_edges = np.histogram(ppts, bins=80, range=(20, 100))
    width = bin_edges[2] - bin_edges[1]
    bin_centres = bin_edges[:-1] + width/2
    plt.bar(bin_centres, freq/run_time, color="C1", alpha=0.7, width=width)
    plt.xlabel(r"PPT /$\mu$s")
    plt.xlim(40,100)
    plt.ylabel("Counts per minute")
    plt.title("Plasma Propation Time all cells")
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/read_red/{}_all_ppt_times.pdf".format(run))


def plot_ppt_all_divided(file_name, run_time):
    run = file_name.split("/")[-1].split("_")[1]
    file = ROOT.TFile(file_name, "READ")
    tree = file.event_tree
    ppts = [[], []]

    for event in tree:
        if len(event.tracker_time_anode) > 50:
            continue
        if len(list(event.calo_time)) == 0:
            continue
        for i in range(len(event.tracker_time_anode)):
            if event.tracker_time_top_cathode[i] != is_none and event.tracker_time_bottom_cathode[i] != is_none and \
                    event.tracker_time_anode[i] != is_none:
                cell_num = int(event.tracker_cell_num[i])
                side, row, layer = cell_id(cell_num)
                if row == 42 or row == 55:
                    index = 1
                elif layer == 0 or layer == 8:
                    index = 1
                else:
                    index = 0
                time = np.amin(list(event.calo_time)) * 1e6
                r5 = event.tracker_time_bottom_cathode[i] * 1e6 - time
                r6 = event.tracker_time_top_cathode[i] * 1e6 - time
                ppts[index].append(r5 + r6)
    file.Close()

    plt.figure(figsize=figsize)
    freq, bin_edges = np.histogram(ppts[0], bins=80, range=(20, 100))
    width = bin_edges[2] - bin_edges[1]
    bin_centres = bin_edges[:-1] + width / 2
    plt.bar(bin_centres, freq / run_time, color="C0", alpha=0.7, width=width, label='Centre Cells')
    freq, bin_edges = np.histogram(ppts[1], bins=80, range=(20, 100))
    width = bin_edges[2] - bin_edges[1]
    bin_centres = bin_edges[:-1] + width / 2
    plt.bar(bin_centres, freq / run_time, color="C1", alpha=0.7, width=width, label='Edge Cells')
    plt.xlabel(r"PPT /$\mu$s")
    plt.xlim(40, 100)
    plt.ylabel("Counts per minute")
    plt.title("Plasma Propation Time all cells")
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/read_red/{}_all_ppt_divided_times.pdf".format(run))


def plot_event_map(file_name, run_time):
    run = file_name.split("/")[-1].split("_")[1]
    file = ROOT.TFile(file_name, "READ")
    tree = file.event_tree
    hits = [None for i in range(2034)]
    for event in tree:
        if len(event.tracker_time_anode) > 50:
            continue
        if len(list(event.calo_time)) == 0:
            continue
        for i in range(len(event.tracker_time_anode)):
            if event.tracker_time_top_cathode[i] != is_none and event.tracker_time_bottom_cathode[i] != is_none and \
                    event.tracker_time_anode[i] != is_none:
                cell_num = int(event.tracker_cell_num[i])
                if hits[cell_num] is None:
                    hits[cell_num] = 1
                else:
                    hits[cell_num] += 1
    file.Close()

    sntracker = sn.tracker('{}_event_map'.format(run), with_palette=True)
    sntracker.draw_cellid = False
    sntracker.draw_cellnum = False
    sntracker.draw_content = False

    for i_cell in range(len(hits)):
        if hits[i_cell] is None:
            continue
        sntracker.setcontent(i_cell, hits[i_cell]/run_time)

    # sntracker.setrange(0, 35)

    sntracker.draw()
    sntracker.save("/Users/williamquinn/Desktop/read_red")


def plot_av_ppt_HV_map():
    files = ['snemo_run-617_red.root', 'snemo_run-635_red.root', 'snemo_run-638_red.root', 'snemo_run-645_red.root',
             'snemo_run-652_red.root', 'snemo_run-660_red.root', 'snemo_run-668_red.root', 'snemo_run-608_red.root']
    HVs = [1600, 1650, 1612, 1625, 1650, 1675, 1600, 1600]
    Areas = [3, 5, 6, 7, 1, 2, 4, 0]
    av_ppt = [None for i in range(2034)]
    hvs = [None for i in range(2034)]

    for index, file in enumerate(files):
        print(f">>> File to process {file}")
        data = parse_root_file("/Users/williamquinn/Desktop/read_red/" + file)
        ppts = [[] for i in range(2034)]
        for event in data:
            cell_data = data[event]["cells"]
            for cell in cell_data:
                ppt = cell_data[cell]["ppts"]
                ts = cell_data[cell]['ts']
                if ppt is None:
                    continue
                elif ts[0] is None:
                    continue
                else:
                    ppts[cell].append(ppt)
        for cell in range(2034):
            if len(ppts[cell]) == 0:
                continue
            elif av_ppt[cell] is not None:
                # This shouldn't happen really
                print("It happened")
                continue
            else:
                av_ppt[cell] = np.average(ppts[cell])
                hvs[cell] = HVs[index]

    sntracker = sn.tracker("av_ppt", with_palette=True)
    sntracker.draw_cellid = False
    sntracker.draw_cellnum = False
    sntracker.draw_content = False
    for cell in range(2034):
        if av_ppt[cell] is None:
            continue
        else:
            sntracker.setcontent(cell, av_ppt[cell])
    sntracker.draw()
    sntracker.save("/Users/williamquinn/Desktop/read_red")
    del sntracker

    sntracker = sn.tracker("full_HV", with_palette=True)
    sntracker.draw_cellid = False
    sntracker.draw_cellnum = False
    sntracker.draw_content = False
    for cell in range(2034):
        if hvs[cell] is None:
            continue
        else:
            sntracker.setcontent(cell, hvs[cell])
    sntracker.draw()
    sntracker.save("/Users/williamquinn/Desktop/read_red")
    del sntracker


def main():
    # args = parse_arguments()
    # input_file = args.i
    # output_directory = args.d

    # plot_full_event_map("/Users/williamquinn/Desktop/read_red", False, output_directory)
    # plot_full_event_map("/Users/williamquinn/Desktop/read_red", True, output_directory)
    # ppts = plot_full_drift_times("/Users/williamquinn/Desktop/read_red", 'full_drift_times', output_directory)
    # plot_drift_time_vs_voltage("/Users/williamquinn/Desktop/read_red", output_directory)
    # plot_ppt_vs_dt("/Users/williamquinn/Desktop/read_red", output_directory)

    # plot_r56_vs_r12("/Users/williamquinn/Desktop/read_red", 'red_669_output.root', output_directory, 10)
    # plot_r56_vs_r12("/Users/williamquinn/Desktop/read_red", 'red_618_output.root', output_directory, 10)
    # plot_r56_vs_r12("/Users/williamquinn/Desktop/read_red", 'red_619_output.root', output_directory, 10)
    # plot_3d_event("/Users/williamquinn/Desktop/read_red/red_619_output.root", ppts)

    '''plot_all_times("/Users/williamquinn/Desktop/read_red/red_612_output.root", 60)
    plot_all_times("/Users/williamquinn/Desktop/read_red/red_619_output.root", 10)
    plot_all_times("/Users/williamquinn/Desktop/read_red/red_625_output.root", 15)
    plot_all_times("/Users/williamquinn/Desktop/read_red/red_626_output.root", 15)
    plot_all_times("/Users/williamquinn/Desktop/read_red/red_627_output.root", 15)
    plot_all_times("/Users/williamquinn/Desktop/read_red/red_659_output.root", 15)
    plot_all_times("/Users/williamquinn/Desktop/read_red/red_660_output.root", 15)
    plot_all_times("/Users/williamquinn/Desktop/read_red/red_661_output.root", 15)'''
    # plot_anode_time_comps("/Users/williamquinn/Desktop/read_red/red_619_output.root", 10)
    # plot_ppt_all("/Users/williamquinn/Desktop/read_red/red_619_output.root", 10)
    # plot_event_map("/Users/williamquinn/Desktop/read_red/red_619_output.root", 10)
    # plot_ppt_all_divided("/Users/williamquinn/Desktop/read_red/red_619_output.root", 10)
    plot_av_ppt_HV_map()

    '''file = open(input_file, "r")
    fl = file.readlines()

    files = []
    for line in fl:
        files.append(line.strip())

    run_times = [60, 10, 10, 10, 15, 3.66667, 30, 5, 5, 5, 15, 10, 10, 10, 15, 10, 10, 10, 10]
    for index, i_file in enumerate(files):
        run = i_file.split("_")[1]

        print(">>> Reading file: {}".format(i_file))
        root_file = ROOT.TFile(input_file.split("filenames.txt")[0] + "/" + i_file.strip(), "READ")
        event_tree = root_file.event_tree

        print(">>> RED Events: {}".format(event_tree.GetEntries()))
        # plot_cathode_events(event_tree, True, run, output_directory, run_times[index])
        # plot_anode_events(event_tree, True, run, output_directory, run_times[index])
        # plot_cathode_events(event_tree, False, run, output_directory, run_times[index])
        # plot_anode_events(event_tree, False, run, output_directory, run_times[index])
        plot_drift_times_all(event_tree, run_times[index], output_directory, run)'''


if __name__ == "__main__":
    main()
