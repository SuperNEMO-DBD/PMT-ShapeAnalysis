from src.PMT_Array import PMT_Array
from src.PMT_Object import PMT_Object
from src.PMT_Waveform import PMT_Waveform
from functions.data_reader_functions import process_xml_file_new
from functions.other_functions import io_parse_arguments
import matplotlib.pyplot as plt


def main():
    args = io_parse_arguments()
    input_file = args.i
    output_file = args.o

    waveforms = process_xml_file_new(input_file)

    pmt_array = PMT_Array([1, 1], "test")

    # Set the cuts from a config file. This also sets up the histograms in the array format
    pmt_array.apply_setting("../config_files/pmt_permeation_config_file.txt")

    pmt_array.set_pmt_templates('templates.root', ['Template_Ch0'])
    pmt_array.get_pmt_object_number(0).set_template_bool(True)

    count = 0

    for i, waveform in enumerate(waveforms):
        pmt_waveform = PMT_Waveform(waveform, pmt_array.get_pmt_object_number(0))

        if pmt_waveform.get_pmt_pulse_times() > 0 and pmt_waveform.get_pmt_pulse_charge() < 100:
            plt.plot(pmt_waveform.get_pmt_waveform())
            count += 1

        if count == 10:
            break


if __name__ == '__main__':
    main()
