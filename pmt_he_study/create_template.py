import sys
sys.path.insert(1, '..')

from src.PMT_Classes import *
from functions.other_functions import *
import time, tqdm
from xml.dom import minidom
import xml.etree.ElementTree as ET


def create_xml_file(input_file: str):

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

                if pmt_waveform.done_sweep and len(pmt_waveform.get_pmt_pulse_times()) > 0:

                    pmt_apulse_times = np.array(pmt_waveform.get_pmt_pulse_times())
                    sweep_start = pmt_waveform.get_pmt_object().get_sweep_range()[0]

                    apulse_event = ET.SubElement(data, 'event')
                    apulse_event.set('ID', str(event_index))
                    apulse_event.set('apulse_num', str(len(pmt_apulse_times)))
                    apulse_waveform = ET.SubElement(apulse_event, 'waveform')
                    apulse_waveform.set('CH', str(channel))
                    apulse_waveform.text = str(trace.firstChild.data)
                    apulse_times = ET.SubElement(apulse_event, 'apulse_times')
                    apulse_times.set('CH', str(channel))
                    apulse_times.text = " ".join([str(i) for i in pmt_apulse_times])
                    apulse_shapes = ET.SubElement(apulse_event, 'apulse_shapes')
                    apulse_shapes.set('CH', str(channel))
                    apulse_shapes.text = " ".join([str(i) for i in pmt_waveform.pmt_waveform_sweep_shape[pmt_apulse_times + sweep_start]])
                    apulse_amps = ET.SubElement(apulse_event, 'apulse_amps')
                    apulse_amps.set('CH', str(channel))
                    apulse_amps.text = " ".join([str(i) for i in pmt_waveform.pmt_waveform_sweep_amp[pmt_apulse_times + sweep_start]])

                    counter += 1
                else:
                    pass
                del pmt_waveform

        if counter % 100 == 0:
            print(counter)
        if counter == 1000:
            break

    _pretty_print(data)
    tree = ET.ElementTree(data)
    tree.write("apulses.xml")


def _pretty_print(current, parent=None, index=-1, depth=0):
    for i, node in enumerate(current):
        _pretty_print(node, current, i, depth + 1)
    if parent is not None:
        if index == 0:
            parent.text = '\n' + ('\t' * depth)
        else:
            parent[index - 1].tail = '\n' + ('\t' * depth)
        if index == len(parent) - 1:
            current.tail = '\n' + ('\t' * (depth - 1))


def main():
    args = io_parse_arguments()
    input_file = args.i
    output_file = args.o

    print(">>> Parsing the data file...")
    processing_start = time.time()
    # parse an xml file by name
    tree = ET.ElementTree()
    tree.parse(input_file)
    root = tree.getroot()
    events = root.findall('event')
    parse_time = time.time() - processing_start
    print(">>> File is good. Parse time: %.3f s" % parse_time)
    print(">>> Number of Events: {}".format(len(events)))

    for event_index, event in tqdm.tqdm(enumerate(events)):
        waveform = event.find('waveform')
        apulse_times = event.find('apulse_times')
        apulse_num = event.attrib['apulse_num']

        print(apulse_num)


if __name__ == '__main__':
    # main()
    create_xml_file('/unix/nemo4/PMT_He_Study_nemo4/data/raw_xml_files/200404/A1400_B1400_t1128.xml')
