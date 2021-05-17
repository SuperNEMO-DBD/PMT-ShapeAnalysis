import sys
sys.path.insert(1, '..')

from src.PMT_Array import PMT_Array
from src.PMT_Object import PMT_Object
from src.PMT_Waveform import PMT_Waveform
from functions.data_reader_functions import process_xml_file_new
from functions.other_functions import io_parse_arguments
import matplotlib.pyplot as plt
from xml.dom import minidom
import time as TIME
import tqdm


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
    processing_start = TIME.time()

    # parse an xml file by name
    xml_file = minidom.parse(input_file)
    events = xml_file.getElementsByTagName('event')
    parse_time = TIME.time() - processing_start

    print(">>> File is good. Parse time: %.3f s" % parse_time)
    print(">>> Number of Events: {}".format(len(events)))

    count = 0

    for event_index, event in tqdm.tqdm(enumerate(events)):

        traces = event.getElementsByTagName('trace')
        for trace_index, trace in enumerate(traces):
            channel = int(trace.attributes['channel'].value)

            if channel == 0:
                pmt_waveform = PMT_Waveform(list(trace.firstChild.data.split(" ")), pmt_array.get_pmt_object_number(0))
                print(pmt_waveform.get_pmt_pulse_peak_position(), pmt_waveform.get_pmt_object().get_pulse_time_threshold(),
                      pmt_waveform.get_pmt_pulse_charge())
                if pmt_waveform.done_sweep:
                    print('done')

                if len(pmt_waveform.get_pmt_pulse_times()) > 0 and abs(pmt_waveform.get_pmt_pulse_charge()) < 100:
                    plt.plot(pmt_waveform.get_pmt_waveform())
                    plt.show()
                    count += 1

                if count == 10:
                    break
                del pmt_waveform


if __name__ == '__main__':
    main()
