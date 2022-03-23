import sys

sys.path.insert(1, '../')
from pmt_he_study.models import *
from ReadRED import sndisplay as sn


def main():
    file = ROOT.TFile("/Users/williamquinn/Desktop/commissioning/run_104.root", "READ")
    tree = file.T
    counter = 0
    x = [i*400/1024 for i in range(1024)]
    plt.figure(figsize=figsize)
    for event in tree:
        if counter == 10:
            break
        om = event.OM_ID + 260
        om_id = om_id_string(om)
        charge = event.charge
        # These seem to be faulty
        # amplitude = event.amplitude
        # baseline = event.baseline

        if om_id == 'M:1.0.0' and counter < 10:
            waveform = list(event.waveform)
            baseline = get_baseline(waveform, 100)
            amplitude = get_amplitude(waveform, baseline)
            if -1 * amplitude > 10:
                waveform = np.array(waveform) - baseline
                plt.plot(x, -1* waveform/amplitude, alpha=0.5, linewidth=1)
                counter += 1
    plt.xlabel("Timestamp /ns")
    # plt.ylabel("ADC Voltage /mV")
    plt.xlim(0, 400)
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/commissioning/waveforms_1_0_0.pdf")


if __name__ == "__main__":
    main()
