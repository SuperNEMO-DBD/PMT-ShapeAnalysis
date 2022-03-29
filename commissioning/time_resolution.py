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
        test_err = [np.sqrt(i)/amplitude for i in test]
        test = [i/amplitude for i in test]
        template_err = [np.sqrt(i/500) for i in template]

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
    graph = ROOT.TGraphErrors(len(shapes), array('d', [i for i in x]), array('d', [i for i in shapes]),
                              array('d', [0 for i in x]), array('d', [i for i in shapes_err]))

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

    return pars, errs, chi


def get_pulse_time_og():
    t = 0
    return t


def read_templates():
    templates = []
    with open("/Users/williamquinn/Desktop/commissioning/templates.csv", "r") as file:
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
        om = event.OM_ID

        if np.sum(counter) == 500*250:
            break

        if counter[om] == n:
            continue

        waveform = list(event.waveform)
        baseline = get_baseline(waveform, 100)
        amplitude = -1 * get_amplitude(waveform, baseline) * adc2mv
        peak = get_peak(waveform)

        if amplitude < 150 or not (25 < peak < 500):
            continue

        temp = waveform[peak - int(25 / tdc2ns):peak + int(75 / tdc2ns)]
        for j in range(len(temp)):
            templates[om][j] += (temp[j] - baseline)/amplitude
        counter[om] += 1

    for i in range(len(templates)):
        if counter[i] > 0:
            templates[i] = np.array(templates[i]) / counter[i]

    with open("/Users/williamquinn/Desktop/commissioning/templates.csv", "w") as out_file:
        for i in range(len(templates)):
            string = ''
            for j in range(len(templates[j])):
                string += str(templates[i][j]) + ","
            string += "\n"
            out_file.write(string)

    file.Close()


def get_event_times(oms, file_name, templates):
    file = ROOT.TFile(file_name, "READ")
    tree = file.T

    events = {
        0: {}
    }

    n_events = tree.GetEntries()
    i_e = 0
    n_plot = 0
    init_time = time.time()
    for event in tree:
        i_e += 1
        if i_e % 10000 == 0:
            temp_time = time.time()
            print(i_e, "/", n_events, temp_time - init_time)
            init_time = temp_time

        om = event.OM_ID
        if templates[om] is None:
            continue
        if om not in oms:
            continue
        event_num = event.event_num
        if event_num in events.keys():
            pass
        else:
            events[event_num] = {}
        waveform = list(event.waveform)
        baseline = get_baseline(waveform, 100)
        amplitude = get_amplitude(waveform, baseline)
        peak = get_peak(waveform)
        tdc = event.tdc * tdc2ns

        if -1 * amplitude < 150 or not (25 < peak < 500):
            continue

        if n_plot < 0:
            pars, errs, chi = get_pulse_time_mf((np.array(waveform) - baseline) / amplitude * -1, templates[om], peak, -1 * amplitude, True, n_plot)
            pulse_time_0 = tdc - 400 + pars[2]
            pulse_time_0_err = errs[2]
            n_plot += 1
        else:
            pars, errs, chi = get_pulse_time_mf((np.array(waveform) - baseline) / amplitude * -1, templates[om], peak, -1 * amplitude, False, n_plot)
            pulse_time_0 = tdc - 400 + pars[2]
            pulse_time_0_err = errs[2]

        pulse_time_1 = (event.fall_cell *tdc2ns)/256.0 - 400 + tdc
        # print(event_num, om, tdc, pulse_time_0, pulse_time_1, pulse_time_1-pulse_time_0)

        events[event_num][om] = [pulse_time_0, pulse_time_0_err, pulse_time_1]
    return events


def store_events(events):
    with open("/Users/williamquinn/Desktop/commissioning/events.csv", "w") as out_file:
        for event in events.keys():
            for om in events[event].keys():
                string = '{},{}'.format(event, om)
                for i in range(len(events[event][om])):
                    string += ',' + str(events[event][om][i])
                string += '\n'
                out_file.write(string)


def read_events():
    events = {}
    with open("/Users/williamquinn/Desktop/commissioning/events.csv", "r") as out_file:
        fl = out_file.readlines()
        for index, line in enumerate(fl):
            line_list = line.split(",")
            event_num = int(line_list[0].strip())
            om = int(line_list[1].strip())
            time1 = float(line_list[2].strip())
            time1_err = float(line_list[3].strip())
            time2 = float(line_list[4].strip())

            if event_num in events.keys():
                events[event_num][om] = [time1, time1_err, time2]
            else:
                events[event_num] = {om: [time1, time1_err, time2]}
    return events


def main():
    '''# create_templates("/Users/williamquinn/Desktop/commissioning/run_430.root")
    templates = read_templates()
    filename = "/Users/williamquinn/Desktop/commissioning/run_430.root"
    events = get_event_times(oms, filename, templates)
    store_events(events)'''
    events = read_events()
    counter = [0, 0, 0]
    for event in events.keys():
        if len(events[event].keys()) > 1:
            if om_num_0 in events[event].keys() and om_num_1 in events[event].keys():
                counter[0] += 1
            elif om_num_0 in events[event].keys() and om_num_2 in events[event].keys():
                counter[1] += 1
            elif om_num_1 in events[event].keys() and om_num_2 in events[event].keys():
                counter[2] += 1
    print(counter)


if __name__ == "__main__":
    main()
