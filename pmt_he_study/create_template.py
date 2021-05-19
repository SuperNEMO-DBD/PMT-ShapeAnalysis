import sys
sys.path.insert(1, '..')

from src.PMT_Classes import *
from functions.other_functions import *
import time, tqdm
from xml.dom import minidom
import xml.etree.ElementTree as ET


def create_xml_file():
    pass


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

    data = ET.Element('data')

    print(">>> Parsing the data file...")
    processing_start = time.time()
    # parse an xml file by name
    xml_file = minidom.parse(input_file)
    events = xml_file.getElementsByTagName('event')
    parse_time = time.time() - processing_start
    print(">>> File is good. Parse time: %.3f s" % parse_time)
    print(">>> Number of Events: {}".format(len(events)))

    counter = 0

    for event_index, event in tqdm.tqdm(enumerate(events)):

        traces = event.getElementsByTagName('trace')
        for trace_index, trace in enumerate(traces):
            channel = int(trace.attributes['channel'].value)

            if channel == 0:
                pmt_waveform = PMT_Waveform(list(trace.firstChild.data.split(" ")), pmt_array.get_pmt_object_number(0))

                if not pmt_waveform.done_sweep:
                    continue

                if len(pmt_waveform.get_pmt_pulse_times()) > 0 and counter < 10000:
                    apulse_event = ET.SubElement(data, 'event')
                    apulse_event.set('ID', str(event_index))
                    apulse_event.set('apulse_num', str(len(pmt_waveform.get_pmt_pulse_times())))
                    apulse_waveform = ET.SubElement(apulse_event, 'waveform')
                    apulse_waveform.set('CH', str(channel))
                    print(str(trace.firstChild.data))
                    apulse_waveform.text(str(trace.firstChild.data))
                    apulse_times = ET.SubElement(apulse_event, 'apulse_times')
                    apulse_times.set('CH', str(channel))
                    apulse_times.text(" ".join([str(i) for i in pmt_waveform.get_pmt_pulse_times()]))

                    counter += 1

    tree = ET.ElementTree()
    tree._setroot(data)
    tree.write('apulses.xml')


if __name__ == '__main__':
    main()
