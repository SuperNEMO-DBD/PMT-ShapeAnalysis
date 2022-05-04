import sys

from xml.dom import minidom
from time import time
import json

import matplotlib.pyplot as plt
import numpy as np
import tqdm

sys.path.insert(1, '../..')
from pmt_he_study.models import *
from ReadRED import sndisplay as sn


def parse_data(input_data_file_name, name):
    print(">>> Parsing the data file...")
    processing_start = time()

    # parse an xml file by name
    xml_file = minidom.parse(input_data_file_name)
    events = xml_file.getElementsByTagName('event')
    parse_time = time() - processing_start

    print(">>> File is good. Parse time: %.3f s" % parse_time)
    print(">>> Number of Events: {}".format(len(events)))

    data = {"baselines": [], "amplitude": [], "charge": [], "peak": []}

    for event in tqdm.tqdm(range(len(events))):
        traces = events[event].getElementsByTagName('trace')
        for trace_index, trace in enumerate(traces):
            waveform = [int(i.strip()) for i in trace.firstChild.data.split(" ")[:-1]]
            baseline = float(np.average(waveform[:150]))
            peak = int(np.argmin(waveform))
            amplitude = float(waveform[peak]) - baseline
            charge = float(-1*(np.sum(waveform[peak - 10: peak+20]) - baseline*30)/50)
            data["baselines"].append(baseline)
            data["amplitude"].append(amplitude)
            data["peak"].append(peak)
            data["charge"].append(charge)
            # waveform_data_list.append(waveform)
    with open("/Users/williamquinn/Desktop/PMT_Project/" + name + "_data.json", "w") as out_file:
        json.dump(data, out_file, indent=4)


def read_data(name):
    with open("/Users/williamquinn/Desktop/PMT_Project/" + name + "_data.json", "r") as out_file:
        data = json.load(out_file)
    return data


def plot_charge(gamma_data, electron_data):
    plt.figure(figsize=figsize)
    freq, bin_edges = np.histogram(gamma_data["charge"], range=(0, 60), bins=240)
    width = bin_edges[2] - bin_edges[1]
    bin_centres = bin_edges[:-1] + width/2
    plt.bar(bin_centres, freq, width=width, alpha=0.5, label='gamma')
    freq, bin_edges = np.histogram(electron_data["charge"], range=(0, 60), bins=240)
    width = bin_edges[2] - bin_edges[1]
    bin_centres = bin_edges[:-1] + width / 2
    plt.bar(bin_centres, freq, width=width, alpha=0.5, label='electron')

    plt.xlabel("Charge /pC")
    plt.legend(loc='best')
    # plt.yscale('log')
    plt.xlim(0, 60)
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/PMT_Project/gamma_charge_comp.pdf")


def main():
    # gamma_file = "/Users/williamquinn/Desktop/PMT_Project/A1000_t1211_gamma.xml"
    # electron_file = "/Users/williamquinn/Desktop/PMT_Project/A1000_t1524.xml"
    # parse_data(gamma_file, "gamma")
    # parse_data(electron_file, "electron")
    gamma_data = read_data("gamma")
    electron_data = read_data("electron")
    plot_charge(gamma_data=gamma_data, electron_data=electron_data)
    print(">>> Got data")


if __name__ == "__main__":
    main()
