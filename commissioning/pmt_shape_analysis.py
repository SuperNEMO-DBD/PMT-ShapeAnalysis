import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(1, '../')
from pmt_he_study.models import *
from ReadRED import sndisplay as sn

tdc2ns = 0.390625
adc2mv = 0.610352


def get_template():
    template = []
    with open("/Users/williamquinn/Desktop/commissioning/template_1_0_1_run_104.csv", "r") as temp_file:
        fl = temp_file.readlines()
        for index, line in enumerate(fl):
            line_list = line.split(",")
            template.append(float(line_list[0].strip()))
    return np.array(template)


def plot_shape_dist(om_id, shapes, run):
    shapes = np.array(shapes)
    print(len(shapes[shapes >= 0.9]), len(shapes[shapes < 0.9]))
    print(np.average(shapes[shapes > 0.9]), np.average(shapes))
    plt.figure(figsize=figsize)
    freq, bin_edges = np.histogram(shapes, range=(0.90, 1), bins=20)
    width = bin_edges[2] - bin_edges[1]
    bin_centres = bin_edges[:-1] + width / 2
    plt.bar(bin_centres, freq, width=width)
    plt.xlabel('MF Shape Index')
    plt.axvline(np.average(shapes), ls='--', color='k',
                label='Mean = {:.4f}'.format(np.average(shapes)))
    plt.xlim(0.9, 1)
    plt.legend(loc='best')
    plt.ylabel("Counts")
    plt.title("PMT Pulse shape index OM - " + om_id)
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/commissioning/{}_shape_dist_{}.pdf".format(run, om_id))


def plot_amp_vs_shape(amps, shapes, run, om_id):
    plt.figure(figsize=figsize)
    plt.plot(amps, shapes, 'C0.', alpha=0.5)
    plt.xlabel("Amplitude /mV")
    plt.ylabel("Shape Index")
    plt.title("PMT Pulse Amplitude vs Shape Index - " + om_id)
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/commissioning/{}_amp_vs_shape_{}.pdf".format(run, om_id))


def plot_amp_vs_fwhm(amps, fwhms, run, om_id):
    plt.figure(figsize=figsize)
    plt.plot(amps, fwhms, 'C1.', alpha=0.5)
    plt.xlabel("Amplitude /mV")
    plt.ylabel("FWHM /ns")
    plt.title("PMT Pulse Amplitude vs Pulse Width - " + om_id)
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/commissioning/{}_amp_vs_fwhm_{}.pdf".format(run, om_id))


def plot_fwhm_vs_shape(fwhms, shapes, run, om_id):
    plt.figure(figsize=figsize)
    plt.plot(fwhms, shapes, 'C3.', alpha=0.5)
    plt.xlabel("FWHM /ns")
    plt.ylabel("Shape Index")
    plt.title("PMT Pulse Width vs Shape Index - " + om_id)
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/commissioning/{}_fwhm_vs_shape_{}.pdf".format(run, om_id))


def plot_fwhm_shape_vs_amp(amps, shapes, fwhms, run, om_id):
    fig, axes = plt.subplots(nrows=1, ncols=2, sharex=True, figsize=figsize)
    fig.suptitle("PMT Pulse Amplitude vs Shape - " + om_id)
    axes[0].plot(amps, shapes, 'C0.', alpha=0.5, markersize=1.5)
    axes[0].set_ylabel("Shape Index")

    axes[1].plot(amps, fwhms, 'C1.', alpha=0.5, markersize=1.5)
    axes[1].yaxis.tick_right()
    axes[1].yaxis.set_label_position("right")
    axes[1].set_ylabel("FWHM /ns")

    fig.supxlabel("Amplitude /mV")
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/commissioning/{}_fwhm_shape_vs_amp_{}.pdf".format(run, om_id))


def plot_av_shapes(avs, stds, ref_om, run):
    sncalo_0 = sn.calorimeter("average_shape_{}".format(run), with_palette=True)
    sncalo_0.draw_content = False
    sncalo_0.draw_omid = False
    for om, val in enumerate(avs):
        if om == ref_om and ref_om is not None:
            sncalo_0.setcontent(om, 1)
        elif val is None or val < 0.9:
            pass
        else:
            sncalo_0.setcontent(om, val)
    sncalo_0.setrange(0.95, 1)
    sncalo_0.draw()
    sncalo_0.save("/Users/williamquinn/Desktop/commissioning")
    del sncalo_0

    sncalo_1 = sn.calorimeter("std_shape_{}".format(run), with_palette=True)
    sncalo_1.draw_content = False
    sncalo_1.draw_omid = False
    for om, val in enumerate(stds):
        if om == ref_om and ref_om is not None:
            sncalo_1.setcontent(om, 1)
        elif val is None:
            pass
        else:
            sncalo_1.setcontent(om, val)
    sncalo_1.draw()
    sncalo_1.save("/Users/williamquinn/Desktop/commissioning")
    del sncalo_1


def plot_waveform(waveform, amplitude, peak, template, om_id, run, shape):
    plt.figure(figsize=figsize)

    start = peak - int(25 / tdc2ns)
    stop = peak + int(175 / tdc2ns)
    x = [i * 400 / 1024 for i in range(1024)][start:stop]
    plt.plot(x, waveform, ".", label='PMT Pulse', markersize=1.5)
    plt.plot(x, amplitude * template/(-1*np.min(template)), label='Template')
    plt.xlabel("Timestamp /ns")
    plt.ylabel('Voltage /mV')
    plt.title("PMT Pulse - {} Shape Index: {:.4f}".format(om_id, shape))
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/commissioning/{}_pulse_vs_template_{}.pdf".format(run, om_id))


def plot_om_type_shape(avs, run):
    new_avs = [[] for i in range(4)]
    for om, av in enumerate(avs):
        om_t, st = om_type(om)
        new_avs[om_t].append(av)

    plt.figure(figsize=figsize)
    for index in range(4):
        freq, bin_edges = np.histogram(new_avs[index], bins=10, range=(0.95, 1))
        width = bin_edges[2] - bin_edges[1]
        bin_centres = bin_edges[:-1] + width/2
        plt.bar(bin_centres, freq, width=width, alpha=0.3)
    plt.ylabel("No. OMs")
    plt.xlabel("Shape Index")
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/commissioning/{}_om_tupe_shape_dist.pdf".format(run))


def main():
    args = io_parse_arguments()
    input_file = args.i
    run = int(input_file.split("/")[-1].split(".")[0].split("run_")[-1])

    template = get_template()

    file = ROOT.TFile(input_file, "READ")
    tree = file.T

    n_counter = [0 for i in range(712)]
    shapes = [[] for i in range(712)]
    mf_amplitudes = [[] for i in range(712)]
    amplitudes = [[] for i in range(712)]
    fwhms = [[] for i in range(712)]

    n_events = tree.GetEntries()
    i_e = 0
    for event in tree:
        i_e += 1
        if i_e % 10000 == 0:
            print(i_e, "/", n_events)

        for index, om in enumerate(list(event.OM_ID)):
            if run == 104:
                om = om + 260

            if om > 712:
                continue

            if n_counter[om] < 1000:
                waveform = list(event.waveform)[index*1024: 1024*(index + 1)]
                baseline = get_baseline(waveform, 100)
                amplitude = get_amplitude(waveform, baseline)
                peak = get_peak(waveform)

                if -1 * amplitude > 50 and 50 < peak < (1024-110):
                    fwhm = get_fwhm_timestamp(waveform, baseline, peak, -1 * amplitude) * tdc2ns
                    x = [i * 400 / 1024 for i in range(1024)][peak - int(25 / tdc2ns):peak + int(175 / tdc2ns)]
                    temp = waveform[peak - int(25 / tdc2ns):peak + int(175 / tdc2ns)]
                    temp = (np.array(temp) - baseline) * adc2mv
                    amplitude = -1* amplitude * adc2mv

                    try:
                        mf_shape, mf_amp = mf_waveform(waveform=temp, template=template)
                    except ValueError:
                        continue
                    if mf_shape < 0.9 and om_id_string(om) == 'M:1.1.1':
                        plot_waveform(temp, amplitude, peak, template, om_id_string(om), run, mf_shape)
                    shapes[om].append(mf_shape)
                    mf_amplitudes[om].append(mf_amp)
                    amplitudes[om].append(amplitude)
                    fwhms[om].append(fwhm)
                    n_counter[om] += 1

    avs = [None for i in range(712)]
    stds = [None for i in range(712)]
    ref_om = None
    for om in range(712):
        om_id = om_id_string(om)

        if len(shapes[om]) == 0:
            av = None
            std = None
        else:
            av = np.average(shapes[om])
            std = np.std(shapes[om])

        avs[om] = av
        stds[om] = std

        if om_id == 'M:1.1.1':
            plot_shape_dist(om_id, shapes[om], run)
            # plot_fwhm_shape_vs_amp(amplitudes[om], shapes[om], fwhms[om], run, om_id)
            # plot_amp_vs_shape(amplitudes[om], shapes[om], run, om_id)
            # plot_amp_vs_fwhm(amplitudes[om], fwhms[om], run, om_id)
            # plot_fwhm_vs_shape(fwhms[om], shapes[om], run, om_id)

        if om_id == 'M:1.0.1':
            ref_om = om

    plot_av_shapes(avs, stds, ref_om, run)
    plot_om_type_shape(avs, run)
    file.Close()


if __name__ == "__main__":
    main()
