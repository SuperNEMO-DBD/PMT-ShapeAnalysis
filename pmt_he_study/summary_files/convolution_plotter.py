import sys

sys.path.insert(1, '..')
import numpy as np
import ROOT
import matplotlib.pyplot as plt
import tqdm
from scipy.optimize import curve_fit
import xml.dom.minidom as minidom
from format_plot import *


def get_baseline(y):
    baseline = 0
    for i in range(100):
        baseline += y[i]
    return baseline / 100


def get_amplitude(y, baseline):
    return np.min(y) - baseline


def get_templates(filename):
    file = ROOT.TFile(filename, "r")
    hist = file.Get("Template_Ch0")
    template = []
    for i in range(hist.GetNbinsX()):
        template.append(hist.GetBinContent(i + 1))
    template = np.array(template)

    template = template/np.sqrt(np.dot(template, template))
    print(len(template))

    '''fig, ax = plt.subplots(figsize=figsize, facecolor='white')
    x = np.array([i for i in range(template.size)])
    plt.plot(x, template, 'C0-', label='PMT Template', linewidth=1)
    plt.xlabel('Timestamp /ns')
    plt.ylabel('Voltage /V')
    plt.xlim(0, template.size)
    # plt.ylim(-1000, 20)
    plt.show()
    plt.close()'''

    return template


def main():
    filename = "/Users/williamquinn/Desktop/PMT_Project/A1400_B1400_t1119.xml"
    file = minidom.parse(filename)
    traces = file.getElementsByTagName('trace')

    template = get_templates("/Users/williamquinn/Desktop/PMT_Project/pmt_short_templates.root")

    for i in tqdm.tqdm(range(int(traces.length))):
        trace = traces[i].firstChild.data.split(" ")[:-1]
        channel = int(traces[i].attributes['channel'].value)

        if channel == 1:
            continue

        waveform = np.array(trace, dtype='float')

        baseline = get_baseline(waveform)
        amplitude = get_amplitude(waveform, baseline)

        if 0 in waveform:
            pass
        else:
            continue

        waveform_r = (waveform - baseline)

        shapes = []
        amps = []
        for i_sweep in range(800, len(waveform) - template.size):
            test = waveform_r[i_sweep: i_sweep + template.size]

            amp = np.dot(test, template)
            shape = amp/np.sqrt(np.dot(test, test))

            '''fig, ax = plt.subplots(figsize=figsize, facecolor='white')
            x = np.array([i + i_sweep for i in range(template.size)])
            ax.plot(x, test, 'C0-', label='Test', linewidth=1)
            ax.set_xlabel('Timestamp /ns')
            ax.set_ylabel('Voltage /V')
            ax.set_xlim(0, template.size)
            #ax.set_ylim(-1000, 20)
            ax2 = ax.twinx()
            ax2.plot(x, template, 'C1-', label='MF Shape', linewidth=1)
            ax2.set_ylabel('Shape Index')
            #ax2.set_ylim(0, 1)
            ax2.set_xlim(0, template.size)
            plt.tight_layout()
            plt.legend(loc='best')
            # fig.savefig("/Users/williamquinn/Desktop/PMT_Project/convolution.pdf")
            plt.show()
            plt.close()'''

            shapes.append(shape)
            amps.append(amp)

        amps = np.array(amps)
        shapes = np.array(shapes)

        if np.max(shapes) > 0.9:

            fig, ax = plt.subplots(figsize=figsize, facecolor='white')
            x = np.array([i for i in range(waveform.size)])
            plt.title("PMT Matched Filter Convolution")
            ax.plot(x, waveform_r, 'C0-', label='PMT Waveform', linewidth=1)
            ax.set_xlabel('Timestamp /ns')
            ax.set_ylabel('Voltage /V')
            #ax.set_xlim(0, waveform.size)
            ax.set_xlim(4000, 4500)
            ax.set_ylim(-55, 5)
            ax2 = ax.twinx()
            ax2.plot(x[800:len(waveform) - template.size], shapes, 'C1-', label='MF Shape', linewidth=1)
            ax2.set_ylabel('Shape Index')
            ax2.set_ylim(-0.2, 1)
            #ax2.set_yscale('log')
            #ax2.set_xlim(0, waveform.size)
            ax2.set_xlim(4000, 4500)
            plt.tight_layout()
            plt.legend(loc='best')
            fig.savefig("/Users/williamquinn/Desktop/PMT_Project/convolution.pdf")
            break


if __name__ == "__main__":
    main()
