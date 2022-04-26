import sys

sys.path.insert(1, '../')
from pmt_he_study.models import *
from ReadRED import sndisplay as sn


adc2mv = 0.610352


def merge_files():

    chain = ROOT.TChain("T")
    for i in [521, 522]:
        file = "/Users/williamquinn/Desktop/commissioning/run_{}_waveform.root".format(i)
        chain.Add(file)

    chain.Merge("/Users/williamquinn/Desktop/commissioning/demon_ap.root")


def plot_event_map(events):
    sncalo = sn.calorimeter("demo_event_map", with_palette=True)
    sncalo.draw_omid = False
    # sncalo.draw_omnum_label()
    sncalo.draw_content = False
    for om, val in enumerate(events):
        if val > 0:
            sncalo.setcontent(om, val)
    sncalo.draw()
    sncalo.save("/Users/williamquinn/Desktop/commissioning")
    del sncalo


def store_waveform(waveform, name, freq):
    with open("/Users/williamquinn/Desktop/commissioning/" + name + ".csv", "w") as file:
        for i, val in enumerate(waveform):
            file.write('{},{}\n'.format(i/freq, val))


def get_temp():
    template = []
    root_file = ROOT.File()


'''def create_templates(filename):
    print(">>> Creating Templates")
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

    file.Close()'''


def main():
    filename = "/Users/williamquinn/Desktop/commissioning/demon_ap.root"
    file = ROOT.TFile(filename, "READ")
    tree = file.T
    events = [0 for i in range(712)]
    i_e = 0
    n_events = tree.GetEntries()
    for event in tree:
        i_e += 1
        if i_e % 10000 == 0:
            print(i_e, "/", n_events)
        for index, om in enumerate(list(event.OM_ID)):
            if events[om] == 2000:
                continue
            waveform = list(event.waveform)[index*1024:1024*(1 + index)]
            baseline = get_baseline(waveform, 20)
            amplitude = get_amplitude(waveform, baseline) * adc2mv * -1
            if amplitude < 300:
                continue
            else:
                # store_waveform((np.array(waveform) - baseline)*adc2mv, "demon_waveform", 0.64)
                events[om] += 1

    file.Close()


if __name__ == "__main__":
    # merge_files()
    main()
