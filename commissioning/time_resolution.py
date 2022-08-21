import sys
import time

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(1, '../')
from pmt_he_study.models import *
from ReadRED import sndisplay as sn
from sklearn.metrics import r2_score


lorentz_string = '[0] * [1] ** 2 / ((x - [2]) ** 2 + [1] ** 2) + [3]'
gaus_string = "[0]*TMath::Gaus(x,[1],[2]) + [3]"

# Define some variables
tdc2ns = 0.390625
adc2mv = 0.610352

# These are the OMs I decided to use
om_num_0, om_num_1, om_num_2 = 227, 214, 226
oms = [om_num_0, om_num_1, om_num_2]
amplitude_cut = 100  #mV


def get_pulse_time_mf(waveform_, template, peak, amplitude, baseline, plot=False, n=0):
    shapes = []
    shapes_err = []
    x = []
    waveform = [(val - baseline)/amplitude for val in waveform_]
    for i in range(peak - 100, peak-20):
        test = waveform[i:i+len(template)]
        # test_err = [np.sqrt(abs(i))/amplitude for i in test]
        # template_err = [np.sqrt(abs(i)/500) for i in template]

        shape = np.dot(test, template)/np.sqrt(np.dot(test, test))
        shapes.append(shape)

        # A = np.sum([template[i]*test[i]*((template_err[i]/template[i])**2 + (test_err[i]/test[i])**2)for i in range(len(template))])/np.dot(template, test)
        # shape_err = shape * np.sqrt(A**2 + (np.sum([i**2 for i in template_err])/(np.dot(template, template)**2)) + (np.sum([i**2 for i in test_err])/(np.dot(test, test)**2)))
        # shapes_err.append(shape_err)

        x.append(i*tdc2ns)
    guess = array('d', [1, float(x[int(len(x)/2)]), 1])
    '''graph = ROOT.TGraphErrors(len(shapes), array('d', [i for i in x]), array('d', [i for i in shapes]),
                              array('d', [0 for i in x]), array('d', [i for i in shapes_err]))'''
    graph = ROOT.TGraph(len(shapes), array('d', [i for i in x]), array('d', [i for i in shapes]))

    # fit = ROOT.TF1("func", "[0]*TMath::Gaus(x,[1],[2])", float(x[0]), float(x[-1]), len(guess))
    fit = ROOT.TF1("func", lorentz_string)
    fit.SetParameter(0, 1)
    fit.SetParameter(1, 10)
    fit.SetParameter(2, float(x[int(len(x)/2)]))
    fit.SetParameter(3, 0.3)

    '''fit.SetParLimits(0, 0, 1)'''
    fit.SetParLimits(1, 0, 100)
    fit.SetParLimits(2, x[0], x[-1])
    #fit.SetParLimits(3, 0, 1)

    graph.Fit("func", "0Q", "", float(x[0]), float(x[-1]))
    fit_result = {
        'pars': [fit.GetParameter(i) for i in range(4)],
        'errs': [fit.GetParError(i) for i in range(4)],
        'chi': fit.GetChisquare(),
        'ndof': fit.GetNDF()
    }
    x = np.array(x)
    y = lorentzian(x, fit_result['pars'][0], fit_result['pars'][2], fit_result['pars'][1]) + fit_result['pars'][3]
    fit_result['R2'] = r2_score(np.array(shapes), y)

    if plot and fit_result['R2'] < 0.99:
        fig1 = plt.figure(figsize=(5, 5), facecolor='w')
        frame1 = fig1.add_axes((.15, .32, .8, .6))
        plt.plot(x, shapes, ".", markersize=3, label='Data')
        plt.plot(x, y, linewidth=1, label='model')
        plt.axvline(np.average(x, weights=shapes), ls='--', color='C2', linewidth=1, label='Data Mean {:.2f} ns'.format(np.average(x, weights=shapes)))
        plt.axvline(fit_result['pars'][2], ls='--', color='C3', linewidth=1, label='Model Mean {:.2f} ns'.format(fit_result['pars'][2]))
        handles, labels = plt.gca().get_legend_handles_labels()
        # patch = patches.Patch(color='white', label=r'$\chi^2_R =$ {:.2f}'.format(fit_result['chi']))
        patch = patches.Patch(color='white', label=r'R$^2=$ {:.2f}'.format(fit_result['R2']))
        handles.extend([patch])
        plt.ylabel("Shape Index")
        plt.xlim(x[0], x[-1])
        plt.ylim(0, 1.1)

        plt.legend(handles=handles, loc='lower left')
        frame2 = fig1.add_axes((.15, .1, .8, .2))
        shapes = np.array(shapes)
        plt.plot(x, (shapes - y) / y, "k.", markersize=3)
        plt.axhline(0, ls='--')
        plt.xlabel('Timestamp /ns')
        plt.ylabel("(model-data)/model")
        plt.xlim(x[0], x[-1])
        plt.tight_layout()

        plt.savefig('/Users/williamquinn/Desktop/commissioning/pulse_time_fit_{}.pdf'.format(n))
        plt.close(fig1)

        fig2 = plt.figure(figsize=figsize)
        plt.plot([i*tdc2ns for i in range(len(waveform))], [val - baseline for val in waveform_])
        plt.xlabel("Timestamp /ns")
        plt.ylabel('Voltage /mV')
        plt.axvline(fit_result['pars'][2], ls='--', color='C3', linewidth=1)
        plt.tight_layout()
        plt.savefig('/Users/williamquinn/Desktop/commissioning/pulse_time_fit_waveform_{}.pdf'.format(n))
        plt.close(fig2)

    del graph
    del fit

    return fit_result


def plot_reflectometry(times, side):

    sncalo = sn.calorimeter("reflectometery_" + side, with_palette=True)
    sncalo.draw_omnum_label()
    for om in times.keys():
        sncalo.setcontent(om, times[om])
    sncalo.draw()
    sncalo.save('/Users/williamquinn/Desktop/commissioning')


def comp_fit_3_gaus(events, case):
    line_styles = ['-', '--']
    plt.figure(figsize=figsize)

    for index, typ in enumerate(events.keys()):
        sel_events = events[typ][case]
        hist = ROOT.TH1D("hist", 'hist', 20, -10, 10)
        for j in sel_events["diff"]:
            hist.Fill(j)

        fit = ROOT.TF1("func", '[0]*TMath::Gaus(x, [1], [2]) + [3]*TMath::Gaus(x, [4], [2]) + [5]*TMath::Gaus(x, [6], [2])')
        fit.SetParameter(0, 100)
        fit.SetParameter(1, -6)
        fit.SetParameter(2, 0.5)
        fit.SetParameter(3, 100)
        fit.SetParameter(4, 0)
        fit.SetParameter(5, 100)
        fit.SetParameter(6, 2.5)
        fit.SetParLimits(0, 0, 1000)
        fit.SetParLimits(1, -10, 10)
        fit.SetParLimits(2, 0, 10)
        fit.SetParLimits(3, 0, 1000)
        fit.SetParLimits(4, -10, 10)
        fit.SetParLimits(5, 0, 1000)
        fit.SetParLimits(6, -10, 10)
        hist.Fit("func", "0Q")

        fit_result = {
            'sig': fit.GetParameter(2),
            'sig_err': fit.GetParError(2),
            'pars': [fit.GetParameter(i) for i in range(7)],
            'errs': [fit.GetParError(i) for i in range(7)],
            'chi': fit.GetChisquare(),
            'ndof': fit.GetNDF()
        }

        x = np.linspace(-10, 10, 100)

        freq, bin_edges = np.histogram(sel_events["diff"], range=(-10, 10), bins=20)
        width = bin_edges[-1] - bin_edges[-2]
        bin_centres = bin_edges[:-1] + width / 2
        plt.bar(bin_centres, freq, width=width, color='C{}'.format(index), alpha=0.3,
                label=typ + r' $\chi^2$/$N_{DoF}$: ' + '{:.2f}/{}'.format(fit_result['chi'], fit_result['ndof']))
        plt.plot(x, np.array(gaussian_noh(x, fit_result['pars'][1], fit_result['pars'][2], fit_result['pars'][0])) + np.array(
            gaussian_noh(x, fit_result['pars'][4], fit_result['pars'][2], fit_result['pars'][3])) + np.array(
            gaussian_noh(x, fit_result['pars'][6], fit_result['pars'][2], fit_result['pars'][5])),
                 'C{}{}'.format(index, line_styles[index]),
                 label=r'$\sigma$ = ({:.2f} ± {:.2f})ns '.format(fit_result['pars'][2], fit_result['errs'][2]))

        del hist
        del fit

    plt.legend(loc='best')
    plt.xlim(-10, 10)
    plt.xlabel('Time Difference /ns')
    plt.tight_layout()
    plt.savefig('/Users/williamquinn/Desktop/commissioning/{}_case_comp.pdf'.format(case))


def comp_fit_2_gaus(events, case):
    line_styles = ['-', '--']
    plt.figure(figsize=figsize)

    for index, typ in enumerate(events.keys()):
        sel_events = events[typ][case]
        hist = ROOT.TH1D("hist", 'hist', 20, -10, 10)
        for j in sel_events["diff"]:
            hist.Fill(j)

        fit = ROOT.TF1("func", '[0]*TMath::Gaus(x, [1], [2]) + [3]*TMath::Gaus(x, [4], [2])')
        fit.SetParameter(0, 100)
        fit.SetParameter(1, -2.5)
        fit.SetParameter(2, 0.5)
        fit.SetParameter(3, 100)
        fit.SetParameter(4, 2.5)
        fit.SetParLimits(0, 0, 1000)
        fit.SetParLimits(1, -10, 10)
        fit.SetParLimits(2, 0, 10)
        fit.SetParLimits(3, 0, 1000)
        fit.SetParLimits(4, -10, 10)
        hist.Fit("func", "0Q")

        fit_result = {
            'sig': fit.GetParameter(2),
            'sig_err': fit.GetParError(2),
            'pars': [fit.GetParameter(i) for i in range(5)],
            'errs': [fit.GetParError(i) for i in range(5)],
            'chi': fit.GetChisquare(),
            'ndof': fit.GetNDF()
        }

        x = np.linspace(-10, 10, 100)

        freq, bin_edges = np.histogram(sel_events["diff"], range=(-10, 10), bins=20)
        width = bin_edges[-1] - bin_edges[-2]
        bin_centres = bin_edges[:-1] + width / 2
        plt.bar(bin_centres, freq, width=width, color='C{}'.format(index), alpha=0.3,
                label=typ + r' $\chi^2$/$N_{DoF}$: ' + '{:.2f}/{}'.format(fit_result['chi'], fit_result['ndof']))
        plt.plot(x, np.array(gaussian_noh(x, fit_result['pars'][1], fit_result['pars'][2], fit_result['pars'][0])) + np.array(gaussian_noh(x, fit_result['pars'][4], fit_result['pars'][2], fit_result['pars'][3])),
                 'C{}{}'.format(index, line_styles[index]),
                 label=r'$\sigma$ = ({:.2f} ± {:.2f})ns '.format(fit_result['pars'][2], fit_result['errs'][2]))

        del hist
        del fit

    plt.legend(loc='best')
    plt.xlim(-10, 10)
    plt.xlabel('Time Difference /ns')
    plt.tight_layout()
    plt.savefig('/Users/williamquinn/Desktop/commissioning/{}_case_comp.pdf'.format(case))


def fit_2_gaus(events, typ, case, plot=False):
    sel_events = events[typ][case]
    hist = ROOT.TH1D("hist", 'hist', 20, -10, 10)
    for i in sel_events["diff"]:
        hist.Fill(i)

    fit = ROOT.TF1("func", '[0]*TMath::Gaus(x, [1], [2]) + [3]*TMath::Gaus(x, [4], [2])')
    fit.SetParameter(0, 100)
    fit.SetParameter(1, -2.5)
    fit.SetParameter(2, 0.5)
    fit.SetParameter(3, 100)
    fit.SetParameter(4, 2.5)
    fit.SetParLimits(0, 0, 1000)
    fit.SetParLimits(1, -10, 10)
    fit.SetParLimits(2, 0, 10)
    fit.SetParLimits(3, 0, 1000)
    fit.SetParLimits(4, -10, 10)
    hist.Fit("func", "0Q")

    fit_result = {
        'sig': fit.GetParameter(2),
        'sig_err': fit.GetParError(2),
        'pars': [fit.GetParameter(i) for i in range(5)],
        'errs': [fit.GetParError(i) for i in range(5)],
        'chi': fit.GetChisquare(),
        'ndof': fit.GetNDF()
    }

    if plot:
        n_sf = get_n_sfs(fit_result['pars'], fit_result['errs'])
        string0 = r'$A =$ {} ± {} $\mu =$ {:.' + str(n_sf[1]) + 'f} ± {:.' + str(n_sf[1]) + 'f} $\sigma =$ {:.' + str(
            n_sf[2]) + 'f} ± {:.' + str(n_sf[2]) + 'f}'
        string1 = r'$A =$ {} ± {} $\mu =$ {:.' + str(n_sf[4]) + 'f} ± {:.' + str(n_sf[4]) + 'f} $\sigma =$ {:.' + str(
            n_sf[2]) + 'f} ± {:.' + str(n_sf[2]) + 'f}'

        x = np.linspace(-10, 10, 100)

        plt.figure(figsize=figsize)
        freq, bin_edges = np.histogram(sel_events["diff"], range=(-10, 10), bins=20)
        width = bin_edges[-1] - bin_edges[-2]
        bin_centres = bin_edges[:-1] + width / 2
        plt.bar(bin_centres, freq, width=width, color='C0', label='data')
        plt.plot(x, np.array(gaussian_noh(x, fit_result['pars'][1], fit_result['pars'][2], fit_result['pars'][0])) + np.array(gaussian_noh(x, fit_result['pars'][4], fit_result['pars'][2], fit_result['pars'][3])),
                 'C1', label='model')
        plt.plot(x, gaussian_noh(x, fit_result['pars'][1], fit_result['pars'][2], fit_result['pars'][0]), 'C2', label=string0.format(int(fit_result['pars'][0]), int(fit_result['pars'][0]),
                                                                                           fit_result['pars'][1], fit_result['pars'][1],
                                                                                           fit_result['pars'][2], fit_result['pars'][2]))
        plt.plot(x, gaussian_noh(x, fit_result['pars'][4], fit_result['pars'][2], fit_result['pars'][3]), 'C3', label=string1.format(int(fit_result['pars'][3]), int(fit_result['errs'][3]),
                                                                                           fit_result['pars'][4], fit_result['errs'][4],
                                                                                           fit_result['pars'][2], fit_result['errs'][2]))
        handles, labels = plt.gca().get_legend_handles_labels()
        patch = patches.Patch(color='white', label=r'$\chi^2 =${:.2f}/{}'.format(fit_result['chi'], fit_result['ndof']))
        handles.extend([patch])
        plt.legend(loc='upper left', handles=handles)
        plt.xlim(-10, 10)
        plt.ylim(0, np.max(freq)*2)
        plt.xlabel('Time Difference /ns')
        plt.tight_layout()
        plt.savefig('/Users/williamquinn/Desktop/commissioning/{}_case_{}.pdf'.format(case, typ))

    del hist
    del fit

    return fit_result


def fit_3_gaus(events, typ, case, plot=False):
    sel_events = events[typ][case]
    hist = ROOT.TH1D("hist", 'hist', 20, -10, 10)
    for i in sel_events["diff"]:
        hist.Fill(i)

    fit = ROOT.TF1("func", '[0]*TMath::Gaus(x, [1], [2]) + [3]*TMath::Gaus(x, [4], [2]) + [5]*TMath::Gaus(x, [6], [2])')
    fit.SetParameter(0, 100)
    fit.SetParameter(1, -6)
    fit.SetParameter(2, 0.5)
    fit.SetParameter(3, 100)
    fit.SetParameter(4, 0)
    fit.SetParameter(5, 100)
    fit.SetParameter(6, 2.5)
    fit.SetParLimits(0, 0, 1000)
    fit.SetParLimits(1, -10, 10)
    fit.SetParLimits(2, 0, 10)
    fit.SetParLimits(3, 0, 1000)
    fit.SetParLimits(4, -10, 10)
    fit.SetParLimits(5, 0, 1000)
    fit.SetParLimits(6, -10, 10)
    hist.Fit("func", "0Q")

    fit_result = {
        'sig': fit.GetParameter(2),
        'sig_err': fit.GetParError(2),
        'pars': [fit.GetParameter(i) for i in range(7)],
        'errs': [fit.GetParError(i) for i in range(7)],
        'chi': fit.GetChisquare(),
        'ndof': fit.GetNDF()
    }

    if plot:
        n_sf = get_n_sfs(fit_result['pars'], fit_result['errs'])
        string0 = r'$A =$ {} ± {} $\mu =$ {:.' + str(n_sf[1]) + 'f} ± {:.' + str(n_sf[1]) + 'f} $\sigma =$ {:.' + str(
            n_sf[2]) + 'f} ± {:.' + str(n_sf[2]) + 'f}'
        string1 = r'$A =$ {} ± {} $\mu =$ {:.' + str(n_sf[4]) + 'f} ± {:.' + str(n_sf[4]) + 'f} $\sigma =$ {:.' + str(
            n_sf[2]) + 'f} ± {:.' + str(n_sf[2]) + 'f}'
        string2 = r'$A =$ {} ± {} $\mu =$ {:.' + str(n_sf[6]) + 'f} ± {:.' + str(n_sf[6]) + 'f} $\sigma =$ {:.' + str(
            n_sf[2]) + 'f} ± {:.' + str(n_sf[2]) + 'f}'

        x = np.linspace(-10, 10, 100)

        plt.figure(figsize=figsize)
        freq, bin_edges = np.histogram(sel_events["diff"], range=(-10, 10), bins=20)
        width = bin_edges[-1] - bin_edges[-2]
        bin_centres = bin_edges[:-1] + width / 2
        plt.bar(bin_centres, freq, width=width, color='C0', label='data')
        plt.plot(x, np.array(gaussian_noh(x, fit_result['pars'][1], fit_result['pars'][2], fit_result['pars'][0])) + np.array(
            gaussian_noh(x, fit_result['pars'][4], fit_result['pars'][2], fit_result['pars'][3])) + np.array(
            gaussian_noh(x, fit_result['pars'][6], fit_result['pars'][2], fit_result['pars'][5])),
                 'C1', label='model')
        plt.plot(x, gaussian_noh(x, fit_result['pars'][1], fit_result['pars'][2], fit_result['pars'][0]), 'C2', label=string0.format(int(fit_result['pars'][0]), int(fit_result['errs'][0]),
                                                                                           fit_result['pars'][1], fit_result['errs'][1],
                                                                                           fit_result['pars'][2], fit_result['errs'][2]))
        plt.plot(x, gaussian_noh(x, fit_result['pars'][4], fit_result['pars'][2], fit_result['pars'][3]), 'C3', label=string1.format(int(fit_result['pars'][3]), int(fit_result['errs'][3]),
                                                                                           fit_result['pars'][4], fit_result['errs'][4],
                                                                                           fit_result['pars'][2], fit_result['errs'][2]))
        plt.plot(x, gaussian_noh(x, fit_result['pars'][6], fit_result['pars'][2], fit_result['pars'][5]), 'C4', label=string2.format(int(fit_result['pars'][5]), int(fit_result['errs'][5]),
                                                                                           fit_result['pars'][6], fit_result['errs'][6],
                                                                                           fit_result['pars'][2], fit_result['errs'][2]))
        handles, labels = plt.gca().get_legend_handles_labels()
        patch = patches.Patch(color='white', label=r'$\chi^2 =${:.2f}/{}'.format(fit_result['chi'], fit_result['ndof']))
        handles.extend([patch])
        plt.legend(loc='upper left', handles=handles)
        plt.xlim(-10, 10)
        plt.ylim(0, np.max(freq) * 7/3)
        plt.xlabel('Time Difference /ns')
        plt.tight_layout()
        plt.savefig('/Users/williamquinn/Desktop/commissioning/{}_case_{}.pdf'.format(case, typ))

    del hist
    del fit

    return fit_result


def read_corrected_times():
    times_corrected = {}

    with open("/Users/williamquinn/Desktop/commissioning/corrected_times_MW_it.txt") as t_file:
        fl = t_file.readlines()
        for index, line in enumerate(fl):
            line_list = line.split(" ")
            side = int(line_list[0])
            col = int(line_list[1])
            row = int(line_list[2])
            time = float(line_list[3].strip())

            om = row + col * 13 + side * 260
            times_corrected[om] = time
    return times_corrected


def read_templates():
    templates = {}
    with open("/Users/williamquinn/Desktop/commissioning/res_templates.csv", "r") as file:
        fl = file.readlines()
        for index, line in enumerate(fl):
            line_list = line.split(",")[:-1]
            temp = []
            for j in range(1, len(line_list)):
                temp.append(float(line_list[j]))
            templates[int(line_list[0])] = np.array(temp)/np.sqrt(np.dot(temp, temp))
    return templates


def plot_templates(templates):
    fig = plt.figure(figsize=figsize)
    for om in templates.keys():
        plt.plot([x*tdc2ns for x in range(len(templates[om]))], -1*templates[om]/min(templates[om]), label=om_id_string(om))
    plt.xlabel("Timestamp /ns")
    plt.legend()
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/commissioning/res_templates.pdf")
    plt.close(fig)


def create_templates(filename, n_average=500):
    print(">>> Create template PMT pulses")
    print(">>> ...")
    file = ROOT.TFile(filename, "READ")
    tree = file.T

    temp_length = 100
    temp_low = 25
    temp_high = 75
    om_count = 3
    counter = [0 for i in range(om_count)]

    # Initialise the templates to be lists of zeros
    templates = [[0.0 for i in range(int(temp_length / tdc2ns))] for j in range(om_count)]

    n_events = tree.GetEntries()
    i_e = 0
    for event in tree:
        i_e += 1
        if i_e % 10000 == 0:
            print(i_e, "/", n_events, np.sum(counter))

        event_oms = list(event.OM_ID)
        for index, om in enumerate(event_oms):

            if counter[oms.index(om)] == n_average:
                continue

            waveform = list(event.waveform)[index*1024:1024*(index + 1)]
            baseline = get_baseline(waveform, 100)
            amplitude = -1 * get_amplitude(waveform, baseline) * adc2mv
            peak = get_peak(waveform)

            if amplitude < 100 or not (25 < peak < 500):
                continue

            temp = waveform[peak - int(temp_low / tdc2ns):peak + int(temp_high / tdc2ns)]
            for j in range(len(temp)):
                templates[oms.index(om)][j] += (temp[j] - baseline)/amplitude
            counter[oms.index(om)] += 1

    for i in range(len(templates)):
        if counter[i] > 0:
            templates[i] = np.array(templates[i]) / counter[i]

    with open("/Users/williamquinn/Desktop/commissioning/res_templates.csv", "w") as out_file:
        for i in range(len(templates)):
            string = f'{oms[i]},'
            for j in range(len(templates[i])):
                string += str(templates[i][j]) + ","
            string += "\n"
            out_file.write(string)

    file.Close()
    print(">>> File written: ", "/Users/williamquinn/Desktop/commissioning/res_templates.csv")


def process_root_file(file_name, templates, corrected_times):
    print(">>> Processing ROOT file:", file_name)
    file = ROOT.TFile(file_name, "READ")
    tree = file.T

    events = {
        "CFD": {
            0: {om_num_0: [], om_num_1: [], "diff": []},
            1: {om_num_0: [], om_num_2: [], "diff": []},
            2: {om_num_1: [], om_num_2: [], "diff": []}
        },
        "MF": {
            0: {om_num_0: [], om_num_1: [], "diff": []},
            1: {om_num_0: [], om_num_2: [], "diff": []},
            2: {om_num_1: [], om_num_2: [], "diff": []}
        }
    }
    n_events = tree.GetEntries()
    n_store = [0, 0, 0]
    i_e = 0
    init_time = time.time()
    for event in tree:
        i_e += 1
        if i_e % 1000 == 0:
            temp_time = time.time()
            print(i_e, "/", n_events, temp_time - init_time, 'events stored:', n_store)

        event_oms = list(event.OM_ID)
        event_waveforms = [list(event.waveform)[index * 1024:1024 * (1 + index)] for index in range(len(event_oms))]
        event_baselines = [get_baseline(waveform, 100) for waveform in event_waveforms]
        event_amplitudes = [-1 * adc2mv * get_amplitude(event_waveforms[index], event_baselines[index]) for index in
                            range(len(event_oms))]
        event_peaks = [get_peak(waveform) for waveform in event_waveforms]
        event_tdcs = list(event.tdc)
        event_fall_times = list(event.fall_time)

        # Redundant but good to keep in just in case (and for running on other files)
        if len(event_oms) != 2 or not all(amp > amplitude_cut for amp in event_amplitudes)\
                or not all(100 < peak < 800 for peak in event_peaks):
            continue

        try:
            case = [x == set(event_oms) for x in
                    [{om_num_0, om_num_1}, {om_num_0, om_num_2}, {om_num_1, om_num_2}]].index(True)
        except ValueError:
            continue

        for index, om in enumerate(event_oms):
            cfd = event_fall_times[index] - 400 + event_tdcs[index] * tdc2ns - corrected_times[om]
            mf_results = get_pulse_time_mf(event_waveforms[index], templates[om],
                                           event_peaks[index], event_amplitudes[index], event_baselines[index],
                                           plot=True, n=i_e)
            mf = mf_results['pars'][2] - 400 + event_tdcs[index] * tdc2ns - corrected_times[om]
            events["CFD"][case][om].append(cfd)
            events["MF"][case][om].append(mf)

        events["CFD"][case]["diff"].append(events["CFD"][case][[*events["CFD"][case].keys()][0]][-1] -
                                           events["CFD"][case][[*events["CFD"][case].keys()][1]][-1])
        events["MF"][case]["diff"].append(events["MF"][case][[*events["MF"][case].keys()][0]][-1] -
                                          events["MF"][case][[*events["MF"][case].keys()][1]][-1])
        n_store[case] += 1

    return events


def store_events(events):
    """
    What the events doc looks like:
    events = {
        "CFD": {
            0: {om_num_0: [], om_num_1: [], "diff": []},
            1: {om_num_0: [], om_num_2: [], "diff": []},
            2: {om_num_1: [], om_num_2: [], "diff": []}
        },
        "MF": {
            0: {om_num_0: [], om_num_1: [], "diff": []},
            1: {om_num_0: [], om_num_2: [], "diff": []},
            2: {om_num_1: [], om_num_2: [], "diff": []}
        }
    }
    """
    with open("/Users/williamquinn/Desktop/commissioning/res_events.csv", "w") as out_file:
        for the_type in events.keys():
            out_file.write(f'the_type,{the_type}\n')
            for case in events[the_type].keys():
                out_file.write(f'case,{case}\n')
                for om in events[the_type][case].keys():
                    out_file.write(f'om,{om}\n')
                    string = 'vec,'
                    for i in range(len(events[the_type][case][om])):
                        string += f'{events[the_type][case][om][i]},'
                    string += '\n'
                    out_file.write(string)
            out_file.write('#\n')


def read_events():
    """
        What the events dict looks like:
        events = {
            "CFD": {
                0: {
                    om_num_0: [],
                     om_num_1: [],
                      "diff": []
                    },
                1: {
                    om_num_0: [],
                     om_num_2: [],
                      "diff": []
                      }, ...
                ...
        }
        """
    events = {}
    with open("/Users/williamquinn/Desktop/commissioning/res_events.csv", "r") as out_file:
        fl = out_file.readlines()
        for index, line in enumerate(fl):
            line_list = line.split(",")
            if '#' in line or line == '\n':
                the_type = None
                case = None
                om = None
                val = None

            if line_list[0] == "the_type":
                the_type = line_list[1].strip()
                events[the_type] = {}
                continue
            elif line_list[0] == 'case':
                case = int(line_list[1].strip())
                events[the_type][case] = {}
                continue
            elif line_list[0] == 'om':
                if line_list[1].strip() == 'diff':
                    om = 'diff'
                else:
                    om = int(line_list[1].strip())
                events[the_type][case][om] = None
                continue
            elif line_list[0] == 'vec':
                val = np.array([float(line_list[i].strip()) for i in range(1, len(line_list[:-1]))], dtype=float)
                events[the_type][case][om] = val
            else:
                pass
    return events


def main():
    corrected_times = read_corrected_times()
    # plot_reflectometry(corrected_times, 'it')

    output_directory = "/Users/williamquinn/Desktop/commissioning/"
    input_filename = "output_run_430_waveform.root"

    # The first time running this will create a file called res_templates.csv
    # create_templates(filename=output_directory+input_filename)
    # Once they are created, no need to keep recreating them
    templates = read_templates()
    # plot_templates(templates)
    filename = output_directory + input_filename

    events = process_root_file(filename, templates, corrected_times)
    store_events(events)
    events = read_events()

    MF_0 = fit_3_gaus(events, 'MF', 0, plot=True)
    CFD_0 = fit_3_gaus(events, 'CFD', 0, plot=True)
    MF_1 = fit_2_gaus(events, 'MF', 1, plot=True)
    CFD_1 = fit_2_gaus(events, 'CFD', 1, plot=True)
    MF_2 = fit_2_gaus(events, 'MF', 2, plot=True)
    CFD_2 = fit_2_gaus(events, 'CFD', 2, plot=True)

    MF = np.array([MF_0['sig'], MF_1['sig'], MF_2['sig']])
    MF_err = np.array([MF_0['sig_err'], MF_1['sig_err'], MF_2['sig_err']])
    CFD = np.array([CFD_0['sig'], CFD_1['sig'], CFD_2['sig']])
    CFD_err = np.array([CFD_0['sig_err'], CFD_1['sig_err'], CFD_2['sig_err']])
    the_matrix = np.array([[1, 1, -1],
                           [1, -1, 1],
                           [-1, 1, 1]])
    the_matrix_0 = np.array([[1, 1, 1],
                             [1, 1, 1],
                             [1, 1, 1]])

    comp_fit_3_gaus(events, 0)
    comp_fit_2_gaus(events, 1)
    comp_fit_2_gaus(events, 2)

    MF_res = np.sqrt(0.5*np.matmul(the_matrix, MF ** 2))
    MF_res_err = np.sqrt(np.dot(MF ** 2, MF_err ** 2))/MF_res
    CFD_res = np.sqrt(0.5 * np.matmul(the_matrix, CFD ** 2))
    CFD_res_err = np.sqrt(np.dot(CFD ** 2, CFD_err ** 2)) / CFD_res

    MF_av = np.average(MF_res, weights=1/MF_res_err**2)
    MF_av_err = np.sqrt(np.average((MF_res - MF_av) ** 2, weights=1 / MF_res_err ** 2)) / np.sqrt(len(MF_res))
    CFD_av = np.average(CFD_res, weights=1 / CFD_res_err**2)
    CFD_av_err = np.sqrt(np.average((CFD_res - CFD_av) ** 2, weights=1 / CFD_res_err**2))/np.sqrt(len(CFD_err))

    print("MF", MF_res, MF_res_err, "av:", MF_av, "err:", MF_av_err)
    print("CFD", CFD_res, CFD_res_err, "av:", CFD_av, "err:", CFD_av_err)


if __name__ == "__main__":
    main()
