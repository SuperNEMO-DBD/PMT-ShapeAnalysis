import sys
sys.path.insert(1, '..')

import ROOT
import time, tqdm
import numpy as np
import matplotlib.pyplot as plt
from functions.other_functions import *
from src.PMT_Classes import *


def main():
    # Handle the file inputs
    args = sncd_parse_arguments()
    input_data_filename = args.i
    output_data_filename = args.o

    # Open file
    file = ROOT.TFile(input_data_filename, "READ")
    file.cd()

    # The tree inside is called "T"
    tree = file.T

    pmt_array = PMT_Array([712, 1], "test")
    # Set the cuts from a config file. This also sets up the histograms in the array format
    pmt_array.apply_setting("../config_files/calo_com_config_file.txt")
    pmt_array.set_pmt_templates('/Users/willquinn/Desktop/snemo_templates.root',
                                ['Template_Ch{}'.format(i) for i in range(712)])

    for i in range(pmt_array.get_pmt_total_number()):
        pmt_array.get_pmt_object_number(i).set_sweep_bool(False)

        pmt_array.get_pmt_object_number(i).set_setting("mf_shape_threshold", 0.8)
        pmt_array.get_pmt_object_number(i).set_setting("mf_amp_threshold", 0.1)

    apulse_nums = []
    # average = [np.array([0.0] * 20) for i in range(712)]
    average = np.array([0.0] * 20)
    apulse_amps = [[] for i in range(712)]
    apulse_shapes = [[] for i in range(712)]

    for event in tqdm.tqdm(tree):
        waveform = np.array(list(event.waveform))
        om = event.OM_ID
        t_id = om_id_string(event.OM_ID)
        apulse_num = event.apulse_num
        apulse_time = event.apulse_times
        main_pulse_time = event.main_pulse_time
        apulse_amp = event.apulse_amplitudes
        apulse_shape = event.apulse_shapes

        baseline = np.average(waveform[:30])

        if int(om) >= 712:
            continue

        if apulse_num > 1:
            continue

        # pmt_waveform = PMT_Waveform(waveform, pmt_array.get_pmt_object_number(om))


        # print(apulse_time)

        for i, t in enumerate(apulse_time):
            if 0 in waveform[t:t+40] or apulse_amp[i] > 60:
                pass
            else:
                # scale = np.min(pmt_waveform.get_pmt_waveform_reduced()[t:t+40]) / np.min(pmt_array.get_pmt_object_number(0).get_template_pmt_pulse())
                # average[om] += pmt_waveform.get_pmt_waveform_reduced()[t:t+20]
                average += waveform[t+3:t + 20+3] - baseline
                '''plt.plot(pmt_waveform.get_pmt_waveform_reduced()[t:t+20])
                plt.title(t_id)
                plt.show(block=False)
                plt.pause(0.01)
                plt.close()'''
                # apulse_shapes[om].append(apulse_shape[i])
                # apulse_amps[om].append(apulse_amp[i])

        # del pmt_waveform

    average = average/np.sqrt(np.sum(average**2))
    average = average - average[0]
    new_file = ROOT.TFile(output_data_filename, "RECREATE")

    for i in range(712):
        name = 'Template_Ch{}'.format(i)
        hist = ROOT.TH1D(name, name, len(average), 0, len(average))
        for j in range(len(average)):
            hist.SetBinContent(1+j, average[j])
        new_file.cd()
        hist.Write()
        del hist
    new_file.Close()


if __name__ == "__main__":
    main()
