import sys
sys.path.insert(1, '..')

from src.PMT_Classes import *
from functions.other_functions import *
import time, tqdm
from xml.dom import minidom
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import ROOT


def create_xml_file(input_file: str):

    pmt_array = PMT_Array([1, 1], "test")
    # Set the cuts from a config file. This also sets up the histograms in the array format
    pmt_array.apply_setting("../config_files/pmt_permeation_config_file.txt")
    pmt_array.set_pmt_templates('/unix/nemo4/PMT_He_Study_nemo4/Templates/new/190621_A1400_B1400_templates.root',
                                ['Template_Ch0'])
    pmt_array.get_pmt_object_number(0).set_sweep_bool(True)

    pmt_array.get_pmt_object_number(0).set_setting("mf_shape_threshold", 0.8)
    pmt_array.get_pmt_object_number(0).set_setting("mf_amp_threshold", 0.1)

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
                    apulse_shapes.text = " ".join([str(i) for i in pmt_waveform.pmt_waveform_sweep_shape[pmt_apulse_times]])
                    apulse_amps = ET.SubElement(apulse_event, 'apulse_amps')
                    apulse_amps.set('CH', str(channel))
                    apulse_amps.text = " ".join([str(i) for i in pmt_waveform.pmt_waveform_sweep_amp[pmt_apulse_times]])

                    counter += 1
                else:
                    pass
                del pmt_waveform

        if counter == 10000:
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

    pmt_array = PMT_Array([2, 1], "test")
    pmt_array.apply_setting("../config_files/pmt_permeation_config_file.txt")
    pmt_array.set_pmt_templates('/Users/willquinn/Desktop/new_template.root',
                                ['Template_Ch0', 'Template_Ch1'])
    pmt_array.get_pmt_object_number(0).set_sweep_bool(True)
    pmt_array.get_pmt_object_number(1).set_sweep_bool(True)

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

    '''plt.plot(pmt_array.get_pmt_object_number(0).get_template_pmt_pulse())
    plt.show()
    plt.plot(pmt_array.get_pmt_object_number(1).get_template_pmt_pulse())
    plt.show()'''

    temp_peak = np.argmin(pmt_array.get_pmt_object_number(0).get_template_pmt_pulse())
    temp_peak_new = np.argmin(pmt_array.get_pmt_object_number(1).get_template_pmt_pulse())

    print(temp_peak, temp_peak_new)

    apulse_nums = []
    average = np.array([0.0]*40)
    tot = 0
    apulse_amps = []
    counter = 0
    counter_ = 0

    for event_index, event in enumerate(events):
        waveform = np.array([int(i) for i in event.find('waveform').text.split(" ")[:-1]])
        pmt_waveform = PMT_Waveform([int(i) for i in event.find('waveform').text.split(" ")[:-1]],
                                    pmt_array.get_pmt_object_number(0))
        pmt_waveform_new = PMT_Waveform([int(i) for i in event.find('waveform').text.split(" ")[:-1]],
                                        pmt_array.get_pmt_object_number(1))

        apulse_time = np.array(pmt_waveform.get_pmt_pulse_times(), dtype='i4')
        apulse_time_new = np.array(pmt_waveform_new.get_pmt_pulse_times(), dtype='i4')
        apulse_amp = pmt_waveform.pmt_waveform_sweep_amp
        apulse_amp_new = pmt_waveform_new.pmt_waveform_sweep_amp
        apulse_shape = pmt_waveform.pmt_waveform_sweep_shape
        apulse_shape_new = pmt_waveform_new.pmt_waveform_sweep_shape
        '''apulse_time = np.array([int(i) for i in event.find('apulse_times').text.split(" ")]) + 800
        apulse_shape = np.array([float(i) for i in event.find('apulse_shapes').text.split(" ")])
        apulse_amp = np.array([float(i) for i in event.find('apulse_amps').text.split(" ")])
        apulse_num = int(event.attrib['apulse_num'])'''


        if len(apulse_time) != len(apulse_time_new):
            if len(apulse_time) > len(apulse_time_new):
                counter += 1
            if len(apulse_time) < len(apulse_time_new):
                counter_ += 1
            #print(len(apulse_time), len(apulse_time_new))
            if len(apulse_time) > len(apulse_time_new):
                for i, t_0 in enumerate(apulse_time):
                    t = 800 + t_0
                    x = [t-40 + i for i in range(120)]
                    fig, ax1 = plt.subplots()
                    ax1.set_xlabel('timestamp (ns)')
                    ax1.set_ylabel('ADC /mV')
                    ax1.tick_params(axis='y')
                    ax1.plot(x, waveform[t - 40:t + 80])
                    ax2 = ax1.twinx()
                    ax2.plot(x, apulse_shape[t_0 - 40 - temp_peak:t_0 + 80 - temp_peak], color='C1', label='short')
                    ax2.plot(x, apulse_shape_new[t_0 - 40 - temp_peak_new:t_0 + 80 - temp_peak_new], color='C2',
                             label='long')
                    ax2.set_ylabel('shape')  # we already handled the x-label with ax1
                    ax2.tick_params(axis='y')
                    ax2.axhline(pmt_array.get_pmt_object_number(0).get_setting("mf_shape_threshold"), 0,
                                8000, ls='--', color='r')
                    '''ax2.plot(t + temp_peak, apulse_shape[t_0], 'xr')
                    ax2.plot(t + temp_peak, apulse_shape_new[t_0], 'xk')'''
                    fig.tight_layout()
                    plt.legend(loc='lower left')
                    plt.show()

            if len(apulse_time_new) > len(apulse_time):
                for i, t_0 in enumerate(apulse_time_new):
                    t = 800 + t_0
                    x = [t - 40 + i for i in range(120)]
                    fig, ax1 = plt.subplots()
                    ax1.set_xlabel('timestamp (ns)')
                    ax1.set_ylabel('ADC /mV')
                    ax1.tick_params(axis='y')
                    ax1.plot(x, waveform[t - 40:t + 80])
                    ax2 = ax1.twinx()
                    ax2.plot(x, apulse_shape[t_0 - 40 - temp_peak:t_0 + 80 - temp_peak], color='C1', label='short')
                    ax2.plot(x, apulse_shape_new[t_0 - 40 - temp_peak_new:t_0 + 80 - temp_peak_new], color='C2', label='long')
                    ax2.set_ylabel('shape')  # we already handled the x-label with ax1
                    ax2.tick_params(axis='y')
                    ax2.axhline(pmt_array.get_pmt_object_number(0).get_setting("mf_shape_threshold"), 0,
                                8000, ls='--', color='r')
                    '''ax2.plot(t + temp_peak_new, apulse_shape[t_0], 'xr')
                    ax2.plot(t + temp_peak_new, apulse_shape_new[t_0 - abs(temp_peak_new - temp_peak)], 'xk')'''
                    fig.tight_layout()
                    plt.legend(loc='lower left')
                    plt.show()
        '''for i,t in enumerate(apulse_time):
            apulse_amps.append(apulse_amp[i])
            if 0 in waveform[t:t+40]:
                continue
            scale = np.min(pmt_waveform.get_pmt_waveform_reduced()[t:t+40]) / np.min(pmt_array.get_pmt_object_number(0).get_template_pmt_pulse())
            average += pmt_waveform.get_pmt_waveform_reduced()[t:t+40]
            tot += 1'''
    print(counter, counter_)

    '''template = pmt_array.get_pmt_object_number(0).get_template_pmt_pulse()

    scale = np.min(average/tot)/np.min(pmt_array.get_pmt_object_number(0).get_template_pmt_pulse())
    plt.plot(average[8:28]/tot, label='new')
    plt.plot(template[8:28] * scale, label='old')
    plt.grid()
    plt.legend(loc='lower right')
    plt.show()

    hist = ROOT.TH1D('Template_Ch0', 'Template_Ch0', len(average[8:28]), 0, len(average[8:28]))
    hist_1 = ROOT.TH1D('Template_Ch1', 'Template_Ch1', len(template), 0, len(template))
    for i in range(len(average[8:28])):
        hist.SetBinContent(i+1, average[8+i])

    for i in range(len(template)):
        hist_1.SetBinContent(i+1, template[i])

    file = ROOT.TFile("~/Desktop/new_template.root", "RECREATE")
    file.cd()
    hist.Write()
    hist_1.Write()

    file.Close()'''


if __name__ == '__main__':
    # main()
    create_xml_file('/unix/nemo4/PMT_He_Study_nemo4/data/raw_xml_files/200404/A1400_B1400_t1128.xml')
