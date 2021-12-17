import ROOT
import numpy as np
import matplotlib.pyplot as plt
import sndisplay as sn
from format_plot import *
from scipy.optimize import curve_fit
from random import randint


is_none = -9999.9
tdc2sec = 1.25e-08


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


def plot_event_map(name: str, hits: list, output_directory):
    sntracker = sn.tracker(name, with_palette=True)
    sntracker.draw_cellid_label()
    sntracker.draw_content_label('{:.2f}')

    for i_cell in range(len(hits)):
        if hits[i_cell] == 0:
            continue
        sntracker.setcontent(i_cell, hits[i_cell])

    # sntracker.setrange(0, 35)

    sntracker.draw()
    sntracker.save(output_directory)


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


def main():
    args = parse_arguments()
    input_file = args.i
    output_directory = args.d

    # plot_full_event_map("/Users/williamquinn/Desktop/read_red", False, output_directory)
    plot_full_event_map("/Users/williamquinn/Desktop/read_red", True, output_directory)
    # plot_full_drift_times("/Users/williamquinn/Desktop/read_red", 'full_drift_times', output_directory)
    # plot_drift_time_vs_voltage("/Users/williamquinn/Desktop/read_red", output_directory)
    # plot_ppt_vs_dt("/Users/williamquinn/Desktop/read_red", output_directory)

    #plot_r56_vs_r12("/Users/williamquinn/Desktop/read_red", 'red_669_output.root', output_directory, 10)
    #plot_r56_vs_r12("/Users/williamquinn/Desktop/read_red", 'red_618_output.root', output_directory, 10)
    #plot_r56_vs_r12("/Users/williamquinn/Desktop/read_red", 'red_619_output.root', output_directory, 10)

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
