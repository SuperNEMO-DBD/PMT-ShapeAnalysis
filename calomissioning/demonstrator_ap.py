import sys

sys.path.insert(1, '../')
from pmt_he_study.models import *
from ReadRED import sndisplay as sn


def main():
    file = ROOT.TFile("/Users/williamquinn/Desktop/PMT_Project/run_521.root", "READ")
    tree = file.T

    aan = [[] for i in range(712)]
    sncalo = sn.calorimeter("aan", True)
    #sncalo.draw_omid_label()
    #sncalo.draw_content_label('{:.3f}')
    sncalo.draw_content = False
    sncalo.draw_omid = False

    for event in tree:
        om = event.OM_ID
        aan[om].append(event.apulse_num)

    new_aan = []
    for om in range(len(aan)):
        if get_pmt_type(om) == 8:
            if len(aan[om]) > 0:
                sncalo.setcontent(om, np.average(aan[om]))
                new_aan.append(np.average(aan[om]))

    sncalo.draw()
    sncalo.save("~/Desktop/PMT_Project")

    del sncalo

    plt.figure(figsize=figsize)
    freq, bin_edges = np.histogram(new_aan, range=(0, 1.2), bins=12)
    width = bin_edges[2] - bin_edges[1]
    bin_centres = bin_edges[:-1] + width/2
    plt.bar(bin_centres, freq, width=width, alpha=0.8, color='C2')
    plt.xlabel("Average After-Pulse Number")
    plt.title("Demonstrator PMT Average After-Pulse per Waveform")
    plt.ylabel("PMT Count")
    plt.xlim(0, 1.2)
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/PMT_Project/aan_demonstrator.pdf")

    file.Close()


if __name__ == "__main__":
    main()
