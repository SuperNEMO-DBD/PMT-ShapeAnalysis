import sys

sys.path.insert(1, '../')
from pmt_he_study.models import *
from ReadRED import sndisplay as sn


tdc2ns = 0.390625
adc2mv = 0.610352


def main():
    file_0 = ROOT.TFile("/Users/williamquinn/Desktop/commissioning/run_104.root", "READ")
    tree_0 = file_0.T
    file_1 = ROOT.TFile("/Users/williamquinn/Desktop/commissioning/run_434.root", "READ")
    tree_1 = file_1.T

    x = [i * 400 / 1024 for i in range(int(200 * 2.56))]
    template_0 = np.array([0.0 for i in range(int(200 * 2.56))])
    template_1 = np.array([0.0 for i in range(int(200 * 2.56))])
    template_2 = np.array([0.0 for i in range(int(200 * 2.56))])
    counter = [0, 0, 0]
    n = 1000
    for event in tree_0:
        for index, om in enumerate(list(event.OM_ID)):
            om = event.OM_ID + 260
            om_id = om_id_string(om)
            if counter[0] == n and counter[1] == n:
                break

            if om_id == 'M:1.0.0' and counter[0] < n:
                waveform = list(event.waveform)[index*1024: 1024*(index + 1)]
                baseline = get_baseline(waveform, 100)
                amplitude = get_amplitude(waveform, baseline)
                peak = get_peak(waveform)
                if -1 * amplitude > 10 and peak < 500:
                    temp = waveform[peak - int(25 * 2.56):peak + int(175 * 2.56)]
                    template_0 += (np.array(temp) - baseline)/amplitude * -1
                    counter[0] += 1

            if om_id == 'M:1.0.1' and counter[1] < n:
                waveform = list(event.waveform)[index*1024: 1024*(index + 1)]
                baseline = get_baseline(waveform, 100)
                amplitude = get_amplitude(waveform, baseline)
                peak = get_peak(waveform)
                if -1 * amplitude > 10 and peak < 500:
                    temp = waveform[peak - int(25 * 2.56):peak + int(175 * 2.56)]
                    template_1 += (np.array(temp) - baseline)/amplitude * -1
                    counter[1] += 1

    for event in tree_1:
        for index, om in enumerate(list(event.OM_ID)):
            om = event.OM_ID
            om_id = om_id_string(om)
            if counter[2] == n:
                break

            if om_id == 'M:1.0.0' and counter[2] < n:
                waveform = list(event.waveform)[index*1024: 1024*(index + 1)]
                baseline = get_baseline(waveform, 100)
                amplitude = get_amplitude(waveform, baseline)
                peak = get_peak(waveform)
                if -1 * amplitude > 10 and peak < 500:
                    temp = waveform[peak - int(25 * 2.56):peak + int(175 * 2.56)]
                    template_2 += (np.array(temp) - baseline) / amplitude * -1
                    counter[2] += 1

    template_0 = template_0 / n
    template_1 = template_1 / n
    template_2 = template_2 / n

    plt.figure(figsize=figsize)
    plt.plot(x, template_0, label='5" PMT Before')
    plt.plot(x, template_2, label='5" PMT After')
    plt.plot(x, template_1, label='8" PMT')
    plt.xlabel("Timestamp /ns")
    plt.xlim(0, 200)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/commissioning/template_comparisons.pdf")


if __name__ == "__main__":
    main()
