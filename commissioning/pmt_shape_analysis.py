import sys

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


def main():
    args = io_parse_arguments()
    input_file = args.i
    run = int(input_file.split("/")[-1].split(".")[0].split("run_")[-1])

    template = get_template()

    file = ROOT.TFile(input_file, "READ")
    tree = file.T

    n_counter = [0 for i in range(712)]
    shapes = [[] for i in range(712)]
    amplitudes = [[] for i in range(712)]

    n_events = tree.GetEntries()
    i_e = 0
    for event in tree:
        i_e += 1
        if i_e % 10000 == 0:
            print(i_e, "/", n_events)

        if run == 104:
            om = event.OM_ID + 260
        else:
            om = event.OM_ID

        if om > 712:
            continue

        if n_counter[om] < 1000:
            waveform = list(event.waveform)
            baseline = get_baseline(waveform, 100)
            amplitude = get_amplitude(waveform, baseline)
            peak = get_peak(waveform)

            if -1 * amplitude > 10:
                x = [i * 400 / 1024 for i in range(1024)][peak - int(25 * 2.56):peak + int(175 * 2.56)]
                temp = waveform[peak - int(25 * 2.56):peak + int(175 * 2.56)]
                temp = (np.array(temp) - baseline) * adc2mv
                amplitude = amplitude * adc2mv

                try:
                    mf_shape, mf_amp = mf_waveform(waveform=temp, template=template)
                except ValueError:
                    continue
                shapes[om].append(mf_shape)
                amplitudes[om].append(mf_amp)
                n_counter[om] += 1

    avs = [None for i in range(712)]
    stds = [None for i in range(712)]
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
            plt.figure(figsize=figsize)
            freq, bin_edges = np.histogram(shapes[om], range=(0.9, 1), bins=20)
            width = bin_edges[2] - bin_edges[1]
            bin_centres = bin_edges[:-1] + width/2
            plt.bar(bin_centres, freq/np.sum(freq), width=width)
            plt.xlabel('MF Shape Index')
            plt.xlim(0.9, 1)
            plt.tight_layout()
            plt.savefig("/Users/williamquinn/Desktop/commissioning/{}_shape_dist_1_1_1.pdf".format(run))

    sncalo_0 = sn.calorimeter("average_shape_{}".format(run), with_palette=True)
    sncalo_0.draw_content = False
    sncalo_0.draw_omid = False
    for om, val in enumerate(avs):
        if val is None or val < 0.9:
            pass
        else:
            sncalo_0.setcontent(om, val)
    sncalo_0.setrange(0.95, 1)
    sncalo_0.draw()
    sncalo_0.save("/Users/williamquinn/Desktop/commissioning")

    sncalo_1 = sn.calorimeter("std_shape_{}".format(run), with_palette=True)
    sncalo_1.draw_content = False
    sncalo_1.draw_omid = False
    for om, val in enumerate(stds):
        if val is None:
            pass
        else:
            sncalo_1.setcontent(om, val)
    sncalo_1.setrange(0.95, 1)
    sncalo_1.draw()
    sncalo_1.save("/Users/williamquinn/Desktop/commissioning")


if __name__ == "__main__":
    main()
