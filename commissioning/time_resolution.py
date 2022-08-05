import sys
import time

sys.path.insert(1, '../')
from pmt_he_study.models import *
from ReadRED import sndisplay as sn

tdc2ns = 0.390625
adc2mv = 0.610352

om_num_0 = 227
om_num_1 = 214
om_num_2 = 226
oms = [om_num_0, om_num_1, om_num_2]


lorentz_string = '[0] * [1] ** 2 / ((x - [2]) ** 2 + [1] ** 2) + [3]'
gaus_string = "[0]*TMath::Gaus(x,[1],[2]) + [3]"


def get_pulse_time_mf(waveform, template, peak, amplitude, plot, n):
    shapes = []
    shapes_err = []
    x = []
    for i in range(peak - 100, peak-20):
        test = waveform[i:i+len(template)]
        test_err = [np.sqrt(abs(i))/amplitude for i in test]
        test = [i/amplitude for i in test]
        template_err = [np.sqrt(abs(i)/500) for i in template]

        shape = np.dot(test, template)/np.sqrt(np.dot(test, test))
        shapes.append(shape)

        A = np.sum([template[i]*test[i]*((template_err[i]/template[i])**2 + (test_err[i]/test[i])**2)for i in range(len(template))])/np.dot(template, test)
        shape_err = shape * np.sqrt(A**2 + (np.sum([i**2 for i in template_err])/(np.dot(template, template)**2)) + (np.sum([i**2 for i in test_err])/(np.dot(test, test)**2)))
        shapes_err.append(shape_err)

        x.append(i * tdc2ns)
    x = np.array(x)
    shapes = np.array(shapes)
    shapes_err = np.array(shapes_err)
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
    pars = []
    errs = []
    chi = fit.GetChisquare() / fit.GetNDF()
    ndof = fit.GetNDF()

    for i in range(4):
        pars.append(fit.GetParameter(i))
        errs.append(fit.GetParError(i))

    if plot:
        fig1 = plt.figure(figsize=(5, 5), facecolor='w')
        frame1 = fig1.add_axes((.15, .32, .8, .6))
        plt.errorbar(x, shapes, yerr=shapes_err, fmt=".", markersize=3,
                     capsize=1, linewidth=1, capthick=1, label='Data')
        y = lorentzian(x, pars[0], pars[2], pars[1]) + pars[3]
        plt.plot(x, lorentzian(x, pars[0], pars[2], pars[1]) + pars[3],
                 linewidth=1, label='model')
        plt.axvline(np.average(x, weights=shapes), ls='--', color='C2', linewidth=1, label='Data Mean {:.2f} ns'.format(np.average(x, weights=shapes)))
        plt.axvline(pars[2], ls='--', color='C3', linewidth=1, label='Model Mean {:.2f} ns'.format(pars[2]))
        handles, labels = plt.gca().get_legend_handles_labels()
        patch = patches.Patch(color='white', label=r'$\chi^2_R =$ {:.2f}'.format(chi))
        handles.extend([patch])
        plt.ylabel("Shape Index")
        plt.xlim(x[0], x[-1])
        plt.ylim(0, 1.1)

        plt.legend(handles=handles, loc='lower left')
        frame2 = fig1.add_axes((.15, .1, .8, .2))
        plt.errorbar(x, (shapes - y) / y, yerr=shapes_err / y, fmt="k.", markersize=3,
                     capsize=1, linewidth=1, capthick=1)
        plt.axhline(0, ls='--')
        plt.xlabel('Timestamp /ns')
        plt.ylabel("(model-data)/model")
        plt.xlim(x[0], x[-1])
        plt.tight_layout()

        plt.savefig('/Users/williamquinn/Desktop/commissioning/pulse_time_fit_{}.pdf'.format(n))

    return pars, errs, chi, ndof


def plot_reflectometry(times, side):

    sncalo = sn.calorimeter("reflectometery_" + side, with_palette=True)
    sncalo.draw_omnum_label()
    for om in times.keys():
        sncalo.setcontent(om, times[om])
    sncalo.draw()
    sncalo.save('/Users/williamquinn/Desktop/commissioning')


def fit_3_gaus(sel_events, case, typ):
    hist = ROOT.TH1D("hist", 'hist', 20, -10, 10)
    for i in sel_events:
        hist.Fill(i)

    fit = ROOT.TF1("func", '[0]*TMath::Gaus(x, [1], [2]) + [3]*TMath::Gaus(x, [4], [2]) + [5]*TMath::Gaus(x, [6], [2])')
    fit.SetParameter(0, 100)
    fit.SetParameter(1, -6)
    fit.SetParameter(2, 0.5)
    fit.SetParameter(3, 100)
    fit.SetParameter(4, 0)
    fit.SetParameter(5, 100)
    fit.SetParameter(7, 2.5)
    fit.SetParLimits(0, 0, 1000)
    fit.SetParLimits(1, -10, 10)
    fit.SetParLimits(2, 0, 10)
    fit.SetParLimits(3, 0, 1000)
    fit.SetParLimits(4, -10, 10)
    fit.SetParLimits(3, 0, 1000)
    fit.SetParLimits(4, -10, 10)
    hist.Fit("func", "0Q")

    pars = []
    errs = []
    chi = fit.GetChisquare() / fit.GetNDF()

    for i in range(7):
        pars.append(fit.GetParameter(i))
        errs.append(fit.GetParError(i))

    n_sf = get_n_sfs(pars, errs)
    string0 = r'$A =$ {} ± {} $\mu =$ {:.' + str(n_sf[1]) + 'f} ± {:.' + str(n_sf[1]) + 'f} $\sigma =$ {:.' + str(
        n_sf[2]) + 'f} ± {:.' + str(n_sf[2]) + 'f}'
    string1 = r'$A =$ {} ± {} $\mu =$ {:.' + str(n_sf[4]) + 'f} ± {:.' + str(n_sf[4]) + 'f} $\sigma =$ {:.' + str(
        n_sf[2]) + 'f} ± {:.' + str(n_sf[2]) + 'f}'
    string1 = r'$A =$ {} ± {} $\mu =$ {:.' + str(n_sf[6]) + 'f} ± {:.' + str(n_sf[6]) + 'f} $\sigma =$ {:.' + str(
        n_sf[2]) + 'f} ± {:.' + str(n_sf[2]) + 'f}'

    x = np.linspace(-10, 10, 100)

    plt.figure(figsize=figsize)
    freq, bin_edges = np.histogram(sel_events, range=(-10, 10), bins=20)
    width = bin_edges[-1] - bin_edges[-2]
    bin_centres = bin_edges[:-1] + width / 2
    plt.bar(bin_centres, freq, width=width, color='C0', label='data')
    plt.plot(x, np.array(gaussian_noh(x, pars[1], pars[2], pars[0])) + np.array(
        gaussian_noh(x, pars[4], pars[2], pars[3])) + np.array(
        gaussian_noh(x, pars[6], pars[2], pars[5])),
             'C1', label='model')
    plt.plot(x, gaussian_noh(x, pars[1], pars[2], pars[0]), 'C2', label=string0.format(int(pars[0]), int(errs[0]),
                                                                                       pars[1], errs[1],
                                                                                       pars[2], errs[2]))
    plt.plot(x, gaussian_noh(x, pars[4], pars[2], pars[3]), 'C3', label=string1.format(int(pars[3]), int(errs[3]),
                                                                                       pars[4], errs[4],
                                                                                       pars[2], errs[2]))
    plt.plot(x, gaussian_noh(x, pars[6], pars[2], pars[5]), 'C4', label=string1.format(int(pars[5]), int(errs[5]),
                                                                                       pars[6], errs[6],
                                                                                       pars[2], errs[2]))
    handles, labels = plt.gca().get_legend_handles_labels()
    patch = patches.Patch(color='white', label=r'$\chi_R =${:.2f}'.format(chi))
    handles.extend([patch])
    plt.legend(loc='upper left', handles=handles)
    plt.xlim(-10, 10)
    plt.ylim(0, np.max(freq) * 7/3)
    plt.xlabel('Time Difference /ns')
    plt.tight_layout()
    plt.savefig('/Users/williamquinn/Desktop/commissioning/{}_case_{}.pdf'.format(case, typ))

    del hist
    del fit

    return pars[2], errs[2]


def comp_fit_3_gaus(events, case):
    line_styles = ['-', '--']
    names = ['MF', 'CFD']
    plt.figure(figsize=figsize)
    for i in range(len(events)):
        hist = ROOT.TH1D("hist", 'hist', 20, -10, 10)
        for j in events[i]:
            hist.Fill(j)

        fit = ROOT.TF1("func", '[0]*TMath::Gaus(x, [1], [2]) + [3]*TMath::Gaus(x, [4], [2]) + [5]*TMath::Gaus(x, [6], [2])')
        fit.SetParameter(0, 100)
        fit.SetParameter(1, -6)
        fit.SetParameter(2, 0.5)
        fit.SetParameter(3, 100)
        fit.SetParameter(4, 0)
        fit.SetParameter(5, 100)
        fit.SetParameter(7, 2.5)
        fit.SetParLimits(0, 0, 1000)
        fit.SetParLimits(1, -10, 10)
        fit.SetParLimits(2, 0, 10)
        fit.SetParLimits(3, 0, 1000)
        fit.SetParLimits(4, -10, 10)
        fit.SetParLimits(3, 0, 1000)
        fit.SetParLimits(4, -10, 10)
        hist.Fit("func", "0Q")

        pars = []
        errs = []
        chi = fit.GetChisquare()
        ndof = fit.GetNDF()

        for j in range(7):
            pars.append(fit.GetParameter(j))
            errs.append(fit.GetParError(j))

        x = np.linspace(-10, 10, 100)

        freq, bin_edges = np.histogram(events[i], range=(-10, 10), bins=20)
        width = bin_edges[-1] - bin_edges[-2]
        bin_centres = bin_edges[:-1] + width / 2
        plt.bar(bin_centres, freq, width=width, color='C{}'.format(i), alpha=0.3, label=names[i] + r' $\chi^2$/$N_{DoF}$: ' + '{:.2f}/{}'.format(chi, ndof))
        plt.plot(x, np.array(gaussian_noh(x, pars[1], pars[2], pars[0])) + np.array(
            gaussian_noh(x, pars[4], pars[2], pars[3])) + np.array(
            gaussian_noh(x, pars[6], pars[2], pars[5])),
                 'C{}{}'.format(i, line_styles[i]),
                 label=r'$\sigma$ = ({:.2f} ± {:.2f})ns '.format(pars[2], errs[2]))

        del hist
        del fit

    plt.legend(loc='best')
    plt.xlim(-10, 10)
    plt.xlabel('Time Difference /ns')
    plt.tight_layout()
    plt.savefig('/Users/williamquinn/Desktop/commissioning/{}_case_comp.pdf'.format(case))


def fit_2_gaus(sel_events, case, typ):
    hist = ROOT.TH1D("hist", 'hist', 20, -10, 10)
    for i in sel_events:
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

    pars = []
    errs = []
    chi = fit.GetChisquare() / fit.GetNDF()

    for i in range(5):
        pars.append(fit.GetParameter(i))
        errs.append(fit.GetParError(i))

    n_sf = get_n_sfs(pars, errs)
    string0 = r'$A =$ {} ± {} $\mu =$ {:.' + str(n_sf[1]) + 'f} ± {:.' + str(n_sf[1]) + 'f} $\sigma =$ {:.' + str(n_sf[2]) + 'f} ± {:.' + str(n_sf[2]) + 'f}'
    string1 = r'$A =$ {} ± {} $\mu =$ {:.' + str(n_sf[4]) + 'f} ± {:.' + str(n_sf[4]) + 'f} $\sigma =$ {:.' + str(n_sf[2]) + 'f} ± {:.' + str(n_sf[2]) + 'f}'

    x = np.linspace(-10, 10, 100)

    plt.figure(figsize=figsize)
    freq, bin_edges = np.histogram(sel_events, range=(-10, 10), bins=20)
    width = bin_edges[-1] - bin_edges[-2]
    bin_centres = bin_edges[:-1] + width / 2
    plt.bar(bin_centres, freq, width=width, color='C0', label='data')
    plt.plot(x, np.array(gaussian_noh(x, pars[1], pars[2], pars[0])) + np.array(gaussian_noh(x, pars[4], pars[2], pars[3])),
             'C1', label='model')
    plt.plot(x, gaussian_noh(x, pars[1], pars[2], pars[0]), 'C2', label=string0.format(int(pars[0]), int(errs[0]),
                                                                                       pars[1], errs[1],
                                                                                       pars[2], errs[2]))
    plt.plot(x, gaussian_noh(x, pars[4], pars[2], pars[3]), 'C3', label=string1.format(int(pars[3]), int(errs[3]),
                                                                                       pars[4], errs[4],
                                                                                       pars[2], errs[2]))
    handles, labels = plt.gca().get_legend_handles_labels()
    patch = patches.Patch(color='white', label=r'$\chi_R =${:.2f}'.format(chi))
    handles.extend([patch])
    plt.legend(loc='upper left', handles=handles)
    plt.xlim(-10, 10)
    plt.ylim(0, np.max(freq)*2)
    plt.xlabel('Time Difference /ns')
    plt.tight_layout()
    plt.savefig('/Users/williamquinn/Desktop/commissioning/{}_case_{}.pdf'.format(case, typ))

    del hist
    del fit

    return pars[2], errs[2]


def comp_fit_2_gaus(events, case):
    line_styles = ['-', '--']
    names = ['MF', 'CFD']
    plt.figure(figsize=figsize)
    for i in range(len(events)):
        hist = ROOT.TH1D("hist", 'hist', 20, -10, 10)
        for j in events[i]:
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

        pars = []
        errs = []

        chi = fit.GetChisquare()
        ndof = fit.GetNDF()

        for j in range(5):
            pars.append(fit.GetParameter(j))
            errs.append(fit.GetParError(j))

        x = np.linspace(-10, 10, 100)

        freq, bin_edges = np.histogram(events[i], range=(-10, 10), bins=20)
        width = bin_edges[-1] - bin_edges[-2]
        bin_centres = bin_edges[:-1] + width / 2
        plt.bar(bin_centres, freq, width=width, color='C{}'.format(i), alpha=0.3, label=names[i] + r' $\chi^2$/$N_{DoF}$: ' + '{:.2f}/{}'.format(chi, ndof))
        plt.plot(x, np.array(gaussian_noh(x, pars[1], pars[2], pars[0])) + np.array(gaussian_noh(x, pars[4], pars[2], pars[3])),
                 'C{}{}'.format(i, line_styles[i]),
                 label=r'$\sigma$ = ({:.2f} ± {:.2f})ns '.format(pars[2], errs[2]))

        del hist
        del fit

    plt.legend(loc='best')
    plt.xlim(-10, 10)
    plt.xlabel('Time Difference /ns')
    plt.tight_layout()
    plt.savefig('/Users/williamquinn/Desktop/commissioning/{}_case_comp.pdf'.format(case))


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


def get_pulse_time_og():
    t = 0
    return t


def read_templates():
    templates = []
    with open("/Users/williamquinn/Desktop/commissioning/res_templates.csv", "r") as file:
        fl = file.readlines()
        for index, line in enumerate(fl):
            line_list = line.split(",")[:-1]
            temp = []
            for j in range(len(line_list)):
                temp.append(float(line_list[j]))
            if temp == [0 for i in range(len(temp))]:
                templates.append(None)
            else:
                templates.append(np.array(temp)/np.sqrt(np.dot(temp, temp)))
    return templates


def create_templates(filename):
    print("Starting")
    file = ROOT.TFile(filename, "READ")
    tree = file.T
    counter = [0 for i in range(712)]
    n = 500
    templates = [[0.0 for i in range(int(100 / tdc2ns))] for j in range(712)]
    n_events = tree.GetEntries()
    i_e = 0
    for event in tree:
        i_e += 1
        if i_e % 10000 == 0:
            print(i_e, "/", n_events, np.sum(counter))
        for index, om in enumerate(list(event.OM_ID)):

            if np.sum(counter) == 500*3:
                break

            if counter[om] == n:
                continue

            waveform = list(event.waveform)[index*1024:1024*(index + 1)]
            baseline = get_baseline(waveform, 100)
            amplitude = -1 * get_amplitude(waveform, baseline) * adc2mv
            peak = get_peak(waveform)

            if amplitude < 100 or not (25 < peak < 500):
                continue

            temp = waveform[peak - int(25 / tdc2ns):peak + int(75 / tdc2ns)]
            for j in range(len(temp)):
                templates[om][j] += (temp[j] - baseline)/amplitude
            counter[om] += 1

    for i in range(len(templates)):
        if counter[i] > 0:
            templates[i] = np.array(templates[i]) / counter[i]

    with open("/Users/williamquinn/Desktop/commissioning/res_templates.csv", "w") as out_file:
        for i in range(len(templates)):
            string = ''
            for j in range(len(templates[j])):
                string += str(templates[i][j]) + ","
            string += "\n"
            out_file.write(string)

    file.Close()


def get_event_times(oms, file_name, templates, corrected_times):
    file = ROOT.TFile(file_name, "READ")
    tree = file.T

    events = {}
    n_events = tree.GetEntries()
    n_store = [0,0,0]
    i_e = 0
    n_plot = 1
    init_time = time.time()
    for event in tree:
        i_e += 1
        if i_e % 1000 == 0:
            temp_time = time.time()
            print(i_e, "/", n_events, temp_time - init_time, 'events stored:', n_store)
            init_time = temp_time
        e_oms = list(event.OM_ID)
        if len(e_oms) != 2:
            continue

        for index, om in enumerate(list(event.OM_ID)):
            if om not in oms:
                continue
            event_num = event.event_num
            waveform = list(event.waveform)[index*1024:1024*(1 + index)]
            baseline = get_baseline(waveform, 100)
            amplitude = -1*adc2mv*get_amplitude(waveform, baseline)
            peak = get_peak(waveform)
            tdc = list(event.tdc)[index]
            # waveforms = [list(event.waveform)[int(i * 1024):int((i + 1) * 1024)] for i in range(len(oms))]
            # baselines = [get_baseline(waveforms[i], 100) for i in range(len(oms))]
            # amplitudes = [-1*adc2mv*get_amplitude(waveforms[i], baselines[i]) for i in range(len(oms))]
            # peaks = [get_peak(waveforms[i]) for i in range(len(oms))]
            # tdcs = [list(event.tdc)[i] for i in range(len(oms))]

            if amplitude < 100 or not (100 < peak < 500):
                store_event = False
            else:
                store_event = True

            '''if abs(tdcs[0] - tdcs[1]) > 10:
                store_event = False'''

            if store_event:
                if event_num not in events.keys():
                    events[event_num] = {}
                if n_plot < 0:
                    pars, errs, chi, ndof = get_pulse_time_mf((np.array(waveform) - baseline) / amplitude,
                                                                  templates[om], peak, amplitude, True, n_plot)
                    n_plot += 1
                else:
                    pars, errs, chi, ndof = get_pulse_time_mf((np.array(waveform) - baseline) / amplitude,
                                                                  templates[om], peak, amplitude, False, n_plot)
                pulse_time_0 = tdc*tdc2ns - 400 + pars[2] - corrected_times[om]
                # pulse_time_0 = tdc - 400 + pars[2]
                pulse_time_0_err = errs[2]

                pulse_time_1 = (list(event.fall_cell)[index] *tdc2ns)/256.0 - 400 + tdc*tdc2ns - corrected_times[om]
                # pulse_time_1 = (event.fall_cell * tdc2ns) / 256.0 - 400 + tdc
                # print(event_num, om, tdc, pulse_time_0, pulse_time_1, pulse_time_1-pulse_time_0)

                events[event_num][om] = [pulse_time_0, pulse_time_1, pars, errs, chi, ndof]
                if om == om_num_0:
                    n_store[0] += 1
                elif om == om_num_1:
                    n_store[1] += 1
                elif om == om_num_2:
                    n_store[2] += 1
    return events


def store_events(events):
    """
    What the events doc looks like:
    events = {
        event_num: int = {
            OM: int = [mf_time: float, cfd_time: float, mf_pars: list, mf_errs: list, chi2: float, ndof: int]
        }
    }
    """
    with open("/Users/williamquinn/Desktop/commissioning/res_events.csv", "w") as out_file:
        for event in events.keys():
            if events[event] is None:
                continue
            for om in events[event].keys():
                string = '{},{}'.format(event, om)
                for i in range(len(events[event][om])):
                    if i in [0, 1, 4, 5]:
                        string += ',' + str(events[event][om][i])
                    else:
                        for j in range(len(events[event][om][i])):
                            string += ',' + str(events[event][om][i][j])
                string += '\n'
                out_file.write(string)


def read_events():
    events = {}
    with open("/Users/williamquinn/Desktop/commissioning/res_events.csv", "r") as out_file:
        fl = out_file.readlines()
        for index, line in enumerate(fl):
            line_list = line.split(",")
            if len(line_list) < 3:
                continue
            event_num = int(line_list[0].strip())
            om = int(line_list[1].strip())
            time1 = float(line_list[2].strip())
            time2 = float(line_list[3].strip())

            if event_num in events.keys():
                events[event_num][om] = [time1, time2]
            else:
                events[event_num] = {om: [time1, time2]}
    return events


def main():
    corrected_times = read_corrected_times()
    # plot_reflectometry(corrected_times, 'it')

    # create_templates("/Users/williamquinn/Desktop/commissioning/run_430_waveform.root")
    templates = read_templates()
    filename = "/Users/williamquinn/Desktop/commissioning/run_430_waveform.root"
    events = get_event_times(oms, filename, templates, corrected_times)
    store_events(events)
    events = read_events()
    counter = [0, 0, 0]
    selected_events = {0: {0: [], 1: []}, 1: {0: [], 1: []}, 2: {0: [], 1: []}}
    for event in events.keys():
        if len(events[event].keys()) > 1:
            if om_num_0 in events[event].keys() and om_num_1 in events[event].keys():
                diff_0 = events[event][om_num_0][0] - events[event][om_num_1][0]
                diff_1 = events[event][om_num_0][1] - events[event][om_num_1][1]
                selected_events[0][0].append(diff_0)
                selected_events[0][1].append(diff_1)
                counter[0] += 1
            elif om_num_0 in events[event].keys() and om_num_2 in events[event].keys():
                diff_0 = events[event][om_num_0][0] - events[event][om_num_2][0]
                diff_1 = events[event][om_num_0][1] - events[event][om_num_2][1]
                selected_events[1][0].append(diff_0)
                selected_events[1][1].append(diff_1)
                counter[1] += 1
            elif om_num_1 in events[event].keys() and om_num_2 in events[event].keys():
                diff_0 = events[event][om_num_1][0] - events[event][om_num_2][0]
                diff_1 = events[event][om_num_1][1] - events[event][om_num_2][1]
                selected_events[2][0].append(diff_0)
                selected_events[2][1].append(diff_1)
                counter[2] += 1
    print(counter)

    ab_mf = fit_3_gaus(selected_events[0][0], 'mf', 'ab')
    ab_og = fit_3_gaus(selected_events[0][1], 'og', 'ab')
    ac_mf = fit_2_gaus(selected_events[1][0], 'mf', 'ac')
    ac_og = fit_2_gaus(selected_events[1][1], 'og', 'ac')
    bc_mf = fit_2_gaus(selected_events[2][0], 'mf', 'bc')
    bc_og = fit_2_gaus(selected_events[2][1], 'og', 'bc')

    comp_fit_3_gaus([selected_events[0][0], selected_events[0][1]], 'ab')
    comp_fit_2_gaus([selected_events[1][0], selected_events[1][1]], 'ac')
    comp_fit_2_gaus([selected_events[2][0], selected_events[2][1]], 'bc')

    err_a_mf = np.sqrt(0.5 * (ab_mf[0]) ** 4 - (bc_mf[0]) ** 4 + (ac_mf[0]) ** 4)
    err_a_og = np.sqrt(0.5 * (ab_og[0]) ** 4 - (bc_og[0]) ** 4 + (ac_og[0]) ** 4)
    err_b_mf = np.sqrt(0.5 * (ab_mf[0]) ** 4 - (ac_mf[0]) ** 4 + (bc_mf[0]) ** 4)
    err_b_og = np.sqrt(0.5 * (ab_og[0]) ** 4 - (ac_og[0]) ** 4 + (bc_og[0]) ** 4)
    err_c_mf = np.sqrt(0.5 * (bc_mf[0]) ** 4 - (ab_mf[0]) ** 4 + (ac_mf[0]) ** 4)
    err_c_og = np.sqrt(0.5 * (bc_og[0]) ** 4 - (ab_og[0]) ** 4 + (ac_og[0]) ** 4)

    print("A:", err_a_og, err_a_mf, om_num_0)
    print("B:", err_b_og, err_b_mf, om_num_1)
    print("C:", err_c_og, err_c_mf, om_num_2)


if __name__ == "__main__":
    main()
