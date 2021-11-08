import sys
import ROOT
sys.path.insert(1, '..')
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from functions.other_functions import io_parse_arguments
import xml.dom.minidom as minidom


def get_charge(waveform, baseline, peak_cell):

    if peak_cell > 1000:
        return 0

    charge = 0.0
    resistance = 50.0
    start = peak_cell-10
    end = peak_cell+20
    for i in range(start, end):
        charge += (waveform[i] - baseline)
    return -1*charge/resistance


def get_peak_cell(waveform):
    peak_cell = 0
    temp = waveform[0]
    for i in range(0, waveform.size):
        if ( waveform[i] < temp ):
            temp = waveform[i]
            peak_cell = i
    return peak_cell


def get_pulse_par(x, y, baseline):
    x = np.array(x)

    pulse_length = len(x)
    amplitude = np.min(y) - baseline;
    position = np.argmin(y)

    start = 0
    stop = 0
    for i in range(position):

        if (y[i] - baseline) < (0.3 * amplitude):
            start = i
        if (y[i] - baseline) < (0.7 * amplitude):
            stop = i

    rise = stop - start
    onset = start

    for i in range(position, pulse_length):
        if (y[i] - baseline) > (amplitude / np.exp(1.0)):
            stop = i;
    fall = stop - position

    return onset, rise, fall


def get_baseline(y):
    baseline = 0
    for i in range(100):
        baseline += y[i]
    return baseline / 100


def get_amplitude(y, baseline):
    return np.min(y) - baseline


def pulse(x, onset, rise, fall, pos):
    y = []
    for i in range(len(x)):
        if i < onset:
            y.append(0)
        else:
            temp = -np.exp(-(x[i] - pos) / rise) + np.exp(-(x[i] - pos) / fall)
            y.append(temp)
    return y


def gaus(x, A, mu, sig):
    return A*np.exp( -0.5*((x-mu)/sig)**2 )


def gaus_fix_sig(x, A, mu):
    return A*np.exp( -0.5*((x-mu)/2.21924508)**2 )


def store_template(filename: str, traces):
    try:
        open(filename, "r")
        print(">>> Template file already exists: {}".format(filename))
        return
    except IOError:
        pass

    average_waveform = []
    average_counter = 0
    for i in range(int(traces.length)):
        trace = traces[i].firstChild.data.split(" ")[:-1]
        channel = int(traces[i].attributes['channel'].value)
        if channel == 1:
            continue
        waveform = np.array(trace, dtype='float')

        baseline = get_baseline(waveform)
        amplitude = get_amplitude(waveform, baseline)

        if amplitude > -100:
            continue

        if 0 in waveform:
            continue
        else:
            if len(average_waveform) == 0:
                average_waveform = np.zeros_like(waveform[:30])
            else:
                pos = np.argmin(waveform - baseline)
                size = np.min(waveform - baseline)
                if pos > 700:
                    continue
                # plt.plot(-1 * (waveform - baseline)[pos - 10:pos + 20] / size)

                average_counter += 1
                average_waveform += -1 * (waveform - baseline)[pos - 10:pos + 20] / size

                if average_counter == 100:
                    break

        if i % 1000 == 0:
            print(i)
    template = average_waveform / average_counter

    with open(filename, "w") as out_file:
        for i in range(len(template)):
            line = '{}\n'.format(template[i])
            out_file.write(line)
    print(">>> Stored template file: {}".format(filename))


def store_charges(filename: str, traces, template):

    try:
        open(filename, "r")
        print(">>> Charge file already exists: {}".format(filename))
        return
    except IOError:
        pass

    out_file = open(filename, "w")
    for i in range(int(traces.length)):
        trace = traces[i].firstChild.data.split(" ")[:-1]
        channel = int(traces[i].attributes['channel'].value)
        if channel == 1:
            continue
        waveform = np.array(trace, dtype='float')

        baseline = get_baseline(waveform)
        amplitude = get_amplitude(waveform, baseline)
        peak = get_peak_cell(waveform)

        if i % 10000 == 0:
            print(i)

        if amplitude > -50:
            continue

        if 0 in waveform:
            pass
        else:
            charge = get_charge(waveform, baseline, peak)
            # charge = np.sum((waveform - baseline)[peak - 10:peak + 20]) / 50
            line = '{}\n'.format(charge)
            out_file.write(line)
            continue

        # pos = np.argmin(waveform - baseline)
        pos = np.where(waveform == 0)[0]

        middle = int(len(pos) / 2)  # = int(2.5) = 2
        pos = pos[middle]  # 3
        size = np.min(waveform - baseline)

        waveform_r = (waveform - baseline)[pos - 10:pos + 20]
        pulse = waveform[pos - 10:pos + 20]

        new_pulse = []
        xi = []
        for j in range(len(pulse)):
            if j > 14:
                break
            if pulse[j] == 0:
                continue
            else:
                new_pulse.append(waveform_r[j])
                xi.append(j)

        popt, pcov = curve_fit(f=gaus_fix_sig, xdata=np.array(xi), ydata=new_pulse, bounds=[[-10000, 0], [0, 40]])

        new_temp = -1 * template * popt[0]
        charge = get_charge(new_temp, 0, get_peak_cell(-1*new_temp))
        # charge = np.sum(new_temp) / 50

        line = '{}\n'.format(charge)
        out_file.write(line)
        # charges.append(charge)

    out_file.close()
    print(">>> Stored charge file: {}".format(filename))


def plot_charge(filename: str, charges):
    fig = plt.figure(figsize=(9, 6), facecolor='white')
    freq, bin_edges = np.histogram(np.array(charges), 100, range=(0, 400))
    width = bin_edges[-1] - bin_edges[-2]
    bin_centres = bin_edges[:-1] + width / 2
    plt.bar(bin_centres, freq, width=width, color="blue")
    plt.xlabel('charge /pC')
    plt.ylabel('counts')
    plt.grid()
    # plt.xlim(200, 300)
    plt.savefig(filename)


def plot_fit(filename: str, charges):
    c = ROOT.TCanvas("bi_spectrum", "Bi Spectrum", 800, 600)
    n_bins = 240
    lower = 230
    higher = 300
    min_charge = 0
    max_charge = 400

    freq, bin_edges = np.histogram(charges, bins=n_bins, range=[min_charge, max_charge])
    bin_width = bin_edges[2] - bin_edges[1]

    bin_centres_0 = (bin_edges + bin_width / 2)[:-1]
    low_bin = int(lower / bin_width)
    high_bin = int(higher / bin_width)
    hist = ROOT.TH1D('', '', n_bins, min_charge, max_charge)

    for i in range(len(charges)):
        # print(charges[i])
        hist.Fill(charges[i])
    hist.GetXaxis().SetRangeUser(lower, higher)

    fit = ROOT.TF1("fit",
                   "[0]*"
                   "(7.08*TMath::Gaus(x,[1],[2]) "
                   " + 1.84*TMath::Gaus(x,[1]*(1 + 72.144/975.651),[2]*1.036) "
                   " + 0.44*TMath::Gaus(x,[1]*(1 + 84.154/975.651),[2]*1.042)) ",
                   lower, higher)

    fit.SetParNames("A", "mu", "sigma")

    fit.SetParLimits(0, 0, 10000)
    fit.SetParLimits(1, lower, higher)
    fit.SetParLimits(2, 0.8, 10)
    fit.SetParameters(319, (higher + lower) / 2, 1.09)

    hist.Fit("fit", "S", "", lower, higher)
    hist.Draw()

    c.Update()
    c.Draw()
    c.SaveAs(filename)

    mu = fit.GetParameter(1)
    mu_err = fit.GetParError(1)
    sig = fit.GetParameter(2)
    sig_err = fit.GetParError(2)
    A = fit.GetParameter(0)
    A_err = fit.GetParError(0)
    chi = fit.GetChisquare() / fit.GetNDF()
    print(mu, sig, A, chi)

    return [[mu, mu_err], [sig, sig_err], [A, A_err], chi]


def get_template(filename: str):
    template = []
    with open(filename, "r") as temp_file:
        fl = temp_file.readlines()
        for i, line in enumerate(fl):
            val = float(line.strip())
            template.append(val)
    return np.array(template)


def get_charges(filename: str):
    charges = []
    with open(filename, "r") as charge_file:
        fl = charge_file.readlines()
        for i, line in enumerate(fl):
            val = float(line.strip())
            charges.append(val)
    return charges


def main():
    args = io_parse_arguments()
    input_file = args.i
    output_path = args.o
    filenames = open(input_file, "r")
    for index, i_filename in enumerate(filenames.readlines()):
        filename = i_filename.strip()

        date = filename.split("/")[-2]
        charge_file = output_path + "/{}_1400V_charges.txt".format(date)
        charge_plot_file = output_path + "/{}_1400V_charge.pdf".format(date)
        charge_fit_file = output_path + "/{}_1400V_charge_fit.pdf".format(date)
        template_file = output_path + "/{}_1400V_templates.txt".format(date)

        print("\n\n>>> Input file: {}".format(filename))
        file = minidom.parse(filename)
        print(">>> Parsed file")
        traces = file.getElementsByTagName('trace')

        store_template(template_file, traces)
        template = get_template(template_file)
        print(">>> Got template")

        store_charges(charge_file, traces, template)
        charges = get_charges(charge_file)
        print(">>> Got Charges")

        plot_charge(charge_plot_file, charges)
        pars = plot_fit(charge_fit_file, charges)

        del file
    print(">>> Finished")


if __name__ == "__main__":
    main()
