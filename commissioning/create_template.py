import sys

sys.path.insert(1, '../')
from pmt_he_study.models import *
from ReadRED import sndisplay as sn


tdc2ns = 0.390625
adc2mv = 0.610352


def main():
    file = ROOT.TFile("/Users/williamquinn/Desktop/commissioning/run_104.root", "READ")
    tree = file.T
    counter = 0
    n = 1000
    x = [i*tdc2ns for i in range(int(200/tdc2ns))]
    template = np.array([0.0 for i in range(int(200/tdc2ns))])
    for event in tree:
        if counter == n:
            break
        om = event.OM_ID + 260
        om_id = om_id_string(om)

        if om_id == 'M:1.0.1' and counter < n:
            waveform = list(event.waveform)
            baseline = get_baseline(waveform, 100)
            amplitude = get_amplitude(waveform, baseline)
            peak = get_peak(waveform)

            temp = waveform[peak - int(25/tdc2ns):peak + int(175/tdc2ns)]
            if -1 * amplitude > 10:
                template += (np.array(temp) - baseline)/amplitude * -1
                counter += 1

    template = template/n
    plt.figure(figsize=figsize)
    plt.plot(x, template)
    plt.xlabel("Timestamp")
    plt.xlim(0, 200)
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/commissioning/template.pdf")

    for event in tree:
        om = event.OM_ID + 260
        om_id = om_id_string(om)

        if om_id == 'M:1.0.1':
            waveform = list(event.waveform)
            baseline = get_baseline(waveform, 100)
            amplitude = get_amplitude(waveform, baseline)
            peak = get_peak(waveform)

            x = [i*400/1024 for i in range(1024)][peak - int(25 * 2.56):peak + int(175 * 2.56)]
            temp = waveform[peak - int(25 * 2.56):peak + int(175 * 2.56)]
            if -1 * amplitude > 10:
                temp = (np.array(temp) - baseline) * adc2mv
                amplitude = amplitude * adc2mv

                mf_shape, mf_amp = mf_waveform(waveform=temp, template=template)
                plt.figure(figsize=figsize)
                plt.plot(x, temp, label='PMT Pulse')
                plt.plot(x, template * amplitude * -1, '--', label='Template')
                plt.xlim(x[0], x[-1])
                plt.xlabel('Timestamp /ns')
                plt.ylabel('Voltage /mV')

                handles, labels = plt.gca().get_legend_handles_labels()
                patch = patches.Patch(color='white', label='MF Shape {:.3f}'.format(mf_shape))
                patch_1 = patches.Patch(color='white', label='MF Amp {:.2f}'.format(mf_amp))
                handles.extend([patch, patch_1])
                plt.legend(handles=handles, loc='best')
                plt.tight_layout()
                plt.savefig("/Users/williamquinn/Desktop/commissioning/temp_vs_waveform.pdf")
                break

    with open("/Users/williamquinn/Desktop/commissioning/template_1_0_1_run_104.csv", "w") as out_file:
        for i, val in enumerate(template):
            out_file.write('{}\n'.format(val))


if __name__ == "__main__":
    main()
