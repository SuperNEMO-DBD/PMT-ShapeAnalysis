import sys

sys.path.insert(1, '../')
from pmt_he_study.models import *
from ReadRED import sndisplay as sn

tdc2ns = 0.390625
adc2mv = 0.610352


def main():
    file = ROOT.TFile("/Users/williamquinn/Desktop/commissioning/run_430.root", "READ")
    tree = file.T
    events = [0 for i in range(712)]
    n_events = tree.GetEntries()
    i_e = 0
    for event in tree:
        i_e += 1
        if i_e % 10000 == 0:
            print(i_e, "/", n_events)
        om = event.OM_ID
        amplitude = event.raw_amplitude

        if (-1*amplitude) * adc2mv > 100:
            events[om] += 1

    sncalo = sn.calorimeter("event_map_430", with_palette=True)
    sncalo.draw_omid = False
    sncalo.draw_omnum_label()
    sncalo.draw_content = True

    for om in range(712):
        if events[om] > 0:
            sncalo.setcontent(om, events[om])

    sncalo.draw()
    sncalo.save("/Users/williamquinn/Desktop/commissioning")


if __name__ == "__main__":
    main()
