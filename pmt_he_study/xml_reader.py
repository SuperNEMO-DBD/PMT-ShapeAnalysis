import sys
sys.path.insert(1, '..')

from src.PMT_Array import PMT_Array
from src.PMT_Object import PMT_Object
from src.PMT_Waveform import PMT_Waveform
from functions.data_reader_functions import process_xml_file_new
from functions.other_functions import io_parse_arguments
import matplotlib.pyplot as plt
from xml.dom import minidom
import time
import tqdm
from scipy.signal import find_peaks
import scipy.fft as fft
import numpy as np


def write_to_file(my_list: list, filename: str):
    with open(filename, 'w') as file:
        for entry in my_list:
            file.write(str(entry) + '\n')


def main():
    args = io_parse_arguments()
    input_file = args.i
    output_file = args.o

    pmt_array = PMT_Array([1, 1], "test")

    # Set the cuts from a config file. This also sets up the histograms in the array format
    pmt_array.apply_setting("../config_files/pmt_permeation_config_file.txt")

    pmt_array.set_pmt_templates('/unix/nemo4/PMT_He_Study_nemo4/Templates/new/190621_A1400_B1400_templates.root',
                                ['Template_Ch0'])
    pmt_array.get_pmt_object_number(0).set_sweep_bool(True)

    print(">>> Parsing the data file...")
    processing_start = time.time()

    # parse an xml file by name
    xml_file = minidom.parse(input_file)
    events = xml_file.getElementsByTagName('event')
    parse_time = time.time() - processing_start

    print(">>> File is good. Parse time: %.3f s" % parse_time)
    print(">>> Number of Events: {}".format(len(events)))

    count = False
    _count = False

    template = pmt_array.get_pmt_object_number(0).get_template_pmt_pulse()

    write_to_file(pmt_array.get_pmt_object_number(0).get_template_pmt_pulse(), "template.txt")

    for event_index, event in tqdm.tqdm(enumerate(events)):

        traces = event.getElementsByTagName('trace')
        for trace_index, trace in enumerate(traces):
            channel = int(trace.attributes['channel'].value)

            if channel == 0:
                pmt_waveform = PMT_Waveform(list(trace.firstChild.data.split(" ")), pmt_array.get_pmt_object_number(0))

                if not pmt_waveform.done_sweep:
                    continue

                signal = pmt_waveform.get_pmt_waveform_reduced()

                NFFT = len(signal)
                f_signal = fft.fft(signal, NFFT)
                f_template = fft.fft(template, NFFT)

                f_cor_sig = f_signal * np.conjugate(f_template)
                cor_sig = fft.ifft(f_cor_sig)

                norm = np.sum(template**2) / len(template)
                cor_sig = cor_sig / norm

                peaks, _ = find_peaks(cor_sig[800:], 5, distance=5)

                print(len(pmt_waveform.get_pmt_pulse_times()), len(peaks))

                if len(pmt_waveform.get_pmt_pulse_times()) != len(peaks):
                    x = [i for i in range(800, pmt_waveform.get_pmt_waveform_length())]
                    plt.plot(x, pmt_waveform.get_pmt_waveform_reduced()[800:])
                    plt.plot(x, pmt_waveform.get_pmt_waveform_reduced()[800:][peaks])
                    plt.plot(x, pmt_waveform.get_pmt_waveform_reduced()[800:][pmt_waveform.get_pmt_pulse_times()])
                    plt.show()

                del pmt_waveform

        if count == _count and count:
            break


if __name__ == '__main__':
    main()
