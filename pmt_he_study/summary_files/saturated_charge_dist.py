import sys

sys.path.insert(1, '..')
import numpy as np
import matplotlib.pyplot as plt
import tqdm
from scipy.optimize import curve_fit
import xml.dom.minidom as minidom
from format_plot import *


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
    return A * np.exp(-0.5 * ((x - mu) / sig) ** 2)


def bi_fun(x, A, mu, sig):
    return A * (gaus(x, 7.08, mu, sig) + gaus(x, 1.84, mu * (1 + 72.144 / 975.651), sig * 1.036) + gaus(x, 0.44, mu * (
                1 + 84.154 / 975.651), sig * 1.042))


def gaus_fix_sig_channel_0(x, A, mu):
    return A*np.exp( -0.5*((x-mu)/sig0)**2 )


def gaus_fix_sig_channel_1(x, A, mu):
    return A*np.exp( -0.5*((x-mu)/sig1)**2 )


def create_average_waveform(traces, filename):
    average_waveform = [[], []]
    average_counter = [0, 0]
    fig = plt.figure(figsize=(9, 6), facecolor='white')
    for i in range(int(traces.length)):
        trace = traces[i].firstChild.data.split(" ")[:-1]
        channel = int(traces[i].attributes['channel'].value)
        waveform = np.array(trace, dtype='float')

        baseline = get_baseline(waveform)
        amplitude = get_amplitude(waveform, baseline)

        if average_counter[channel] == 100:
            continue

        if amplitude > -100:
            continue

        if 0 in waveform:
            continue
        else:
            if len(average_waveform[channel]) == 0:
                average_waveform[channel] = np.zeros_like(waveform[:30])
            else:
                pos = np.argmin(waveform - baseline)
                size = np.min(waveform - baseline)
                if pos > 700:
                    continue
                plt.plot(-1 * (waveform - baseline)[pos - 10:pos + 20] / size)

                average_counter[channel] += 1
                average_waveform[channel] += -1 * (waveform - baseline)[pos - 10:pos + 20] / size

                if average_counter[0] == 100 and average_counter[1] == 100:
                    break

        if i % 1000 == 0:
            print(i)

    plt.grid()
    plt.xlabel('timestamp /ns')
    plt.savefig('/Users/williamquinn/Desktop/pulses.pdf')

    out_file = open(filename, "w")
    for i_om in range(len(average_waveform)):
        out_file.write('channel,{}\n'.format(i_om))
        for j in range(len(average_waveform[i_om])):
            out_file.write('{},{}\n'.format(j, average_waveform[i_om][j]/average_counter[i_om]))
    out_file.close()


def get_average_waveform(filename):
    average_waveform = [[], []]
    in_file = open(filename, "r")
    fl = in_file.readlines()
    channel = 0
    for index, line in enumerate(fl):
        line_list = line.split(',')
        if line_list[0].strip() == 'channel':
            channel = int(line_list[-1].strip())
            continue
        average_waveform[channel].append(float(line.split(",")[-1].strip()))
    average_waveform[0] = np.array(average_waveform[0])
    average_waveform[1] = np.array(average_waveform[1])

    return average_waveform


def get_pulse_width(average_waveforms):
    sigmas = []
    for i_om in range(2):
        template = average_waveforms[i_om]

        fig = plt.figure(figsize=figsize, facecolor='white')
        plt.plot(template, ".", label='template')
        x = np.array([i for i in range(len(template))])
        popt, pcov = curve_fit(f=gaus, xdata=x[:-16], ydata=template[:-16],
                               bounds=[[-1000, 0, 0], [0, 1000, 1000]])
        plt.plot(x[:-16], template[:-16], ".", label='points used')
        A = popt[0]
        mu = popt[1]
        sig = popt[2]
        plt.plot(x, gaus(x, *popt))
        plt.xlabel('timestamp /ns')
        plt.title(r'Model: $A=${:.2f} $\mu=${:.2f} $\sigma=${:.2f}'.format(A, mu, sig))
        plt.legend()
        # plt.yscale('log')
        plt.savefig('/Users/williamquinn/Desktop/PMT_Project/pulse_fit_{}.pdf'.format(i_om))
        sigmas.append(sig)
    return sigmas


def write_charges(traces, filename, templates):

    charges = [[], []]
    sat_charges = [[], []]
    ap_charges = [[], []]
    charge_file = open(filename, "w")
    charge_file.write('channel, charge, sat_charge, ap_charge\n')
    for i in tqdm.tqdm(range(int(traces.length))):
        trace = traces[i].firstChild.data.split(" ")[:-1]
        channel = int(traces[i].attributes['channel'].value)

        waveform = np.array(trace, dtype='float')

        baseline = get_baseline(waveform)
        amplitude = get_amplitude(waveform, baseline)

        if amplitude > -50:
            continue

        if 0 in waveform:
            pass
        else:
            charge = np.sum((waveform - baseline)[pos - 10:pos + 20]) / 50
            ap_charge = np.sum((waveform - baseline)[800:]) / 50

            charges[channel].append(charge)
            sat_charges[channel].append(charge)
            ap_charges[channel].append(ap_charge)
            charge_file.write('{},{},{},{}\n'.format(channel, charge, charge, ap_charge))
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

        if channel == 0:
            popt, pcov = curve_fit(f=gaus_fix_sig_channel_0, xdata=np.array(xi), ydata=new_pulse,
                                   bounds=[[-10000, 0], [0, 40]])
        else:
            popt, pcov = curve_fit(f=gaus_fix_sig_channel_1, xdata=np.array(xi), ydata=new_pulse,
                                   bounds=[[-10000, 0], [0, 40]])

        new_temp = -1 * templates[channel] * popt[0]
        sat_charge = np.sum(new_temp) / 50
        charge = np.sum((waveform - baseline)[pos - 10:pos + 20]) / 50
        ap_charge = np.sum((waveform - baseline)[800:]) / 50

        sat_charges[channel].append(sat_charge)
        charges[channel].append(charges)
        ap_charges[channel].append(ap_charge)
        charge_file.write('{},{},{},{}\n'.format(channel, charge, sat_charge, ap_charge))
    charge_file.close()


def get_charges(filename):
    charges = [[], []]
    sat_charges = [[], []]
    ap_charges = [[], []]

    in_file = open(filename, 'r')
    fl = in_file.readlines()
    for index, line in enumerate(fl):
        if index == 0:
            continue
        line_list = line.split(",")
        channel = int(line_list[0])
        charge = float(line_list[1])
        sat_charge = float(line_list[2])
        ap_charge = float(line_list[3])

        charges[channel].append(-1*charge)
        sat_charges[channel].append(-1*sat_charge)
        ap_charges[channel].append(-1*ap_charge)
    return charges, sat_charges, ap_charges


def plot_charge_comp(charges, sat_charges):
    fig = plt.figure(figsize=figsize, facecolor='white')

    ch_charges = np.array(charges)
    ch_sat_charges = np.array(sat_charges)

    lower, upper = 0, 300
    n_bins = 150

    freq, bin_edges = np.histogram(ch_charges, n_bins, range=(lower, upper))
    width = bin_edges[-1] - bin_edges[-2]
    bin_centres = bin_edges[:-1] + width / 2

    sat_freq, sat_bin_edges = np.histogram(ch_sat_charges, n_bins, range=(lower, upper))
    sat_width = sat_bin_edges[-1] - sat_bin_edges[-2]
    sat_bin_centres = sat_bin_edges[:-1] + sat_width / 2

    num = np.sum(ch_charges)
    plt.bar(bin_centres, freq / num, width=width, alpha=0.3, color='C0')
    plt.plot(bin_centres, freq / num, "C0o", markersize=2, label='Charge')
    plt.bar(sat_bin_centres, sat_freq / num, width=width, alpha=0.3, color='C1')
    plt.plot(sat_bin_centres, sat_freq / num, "C1s", markersize=2, label='Reconstruction')

    plt.xlabel('Charge /pC')
    plt.ylabel('Normalised Counts')
    plt.title("1400V PMT Pulse Charge Reconstruction Comparison")
    plt.xlim(lower, upper)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/PMT_Project/charge_comp.pdf")


def plot_tot(traces, templates):
    for i in range(int(traces.length)):
        trace = traces[i].firstChild.data.split(" ")[:-1]
        channel = int(traces[i].attributes['channel'].value)

        waveform = np.array(trace, dtype='float')
        if channel == 1:
            continue

        baseline = get_baseline(waveform)
        amplitude = get_amplitude(waveform, baseline)

        if amplitude > -50:
            continue

        if 0 in waveform:
            pass
        else:
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

        if channel == 0:
            popt, pcov = curve_fit(f=gaus_fix_sig_channel_0, xdata=np.array(xi), ydata=new_pulse,
                                   bounds=[[-10000, 0], [0, 40]])
        else:
            popt, pcov = curve_fit(f=gaus_fix_sig_channel_1, xdata=np.array(xi), ydata=new_pulse,
                                   bounds=[[-10000, 0], [0, 40]])

        new_temp = -1 * templates[channel] * popt[0]

        fig = plt.figure(figsize=figsize, facecolor='white')
        plt.plot([k for k in range(len(waveform_r))], waveform_r, '.', label='PMT Pulse')
        plt.xlabel('Timestamp /ns')
        plt.ylabel('Voltage /mV')
        plt.plot([k for k in range(len(waveform_r))], new_temp, '.', label='Scaled Template')
        plt.plot(np.linspace(0, 15, 100), gaus_fix_sig_channel_0(np.linspace(0, 15, 100), *popt),
                 label='Model')
        plt.plot(xi, new_pulse, '.', label='Model Points')
        plt.legend(loc='best')
        plt.title('Time-Over-Threshold Pulse Reconstruction')
        plt.tight_layout()
        plt.savefig("/Users/williamquinn/Desktop/PMT_Project/tot_plot.pdf")

        break


def main():
    filename = "/Users/williamquinn/Desktop/PMT_Project/A1400_B1400_t1119.xml"
    file = minidom.parse(filename)
    traces = file.getElementsByTagName('trace')

    av_pulse_file = "/Users/williamquinn/Desktop/PMT_Project/average_pulses.csv"
    charge_file = "/Users/williamquinn/Desktop/PMT_Project/charges.csv"

    #create_average_waveform(traces, av_pulse_file)
    average_waveforms = get_average_waveform(av_pulse_file)

    sigmas = get_pulse_width(average_waveforms)
    global sig0
    global sig1
    sig0, sig1 = sigmas[0], sigmas[1]

    plot_tot(traces, average_waveforms)

    #write_charges(traces, charge_file, average_waveforms)
    #charges, sat_charges, ap_charges = get_charges(charge_file)

    #plot_charge_comp(charges[0], sat_charges[0])


if __name__ == "__main__":
    main()
