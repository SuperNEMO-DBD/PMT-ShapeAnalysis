import ROOT
import numpy as np
from scipy.signal import find_peaks


class PMT_Object:

    def __init__(self, pmt_id: str, data_id: str):
        self.pmt_id = pmt_id
        self.data_id = data_id

        self.number_of_events = 0
        self.sweep_bool = False
        self.template_bool = False

        # Default settings
        charge_cut = 6  # pC
        charge_range = [0, 100]
        nbins = 100
        amp_range = [0, 100]
        mf_shape_range = [-1, 1]
        mf_amp_range = [0, 100]
        sweep_range = [0, 500]
        pulse_time_threshold = 100
        apulse_region = 500
        resistance = 50 # Ohms
        mf_shape_threshold = 0.9
        mf_amp_threshold = 25
        baseline = 1000
        waveform_length = 8000
        trigger_point = 100
        integration = [0.3, 0.3]
        he_region = [1400, 2000]

        self.setting_dict = {
            "charge_cut"            : charge_cut,
            "charge_range"          : charge_range,
            "nbins"                 : nbins,
            "amp_range"             : amp_range,
            "mf_shape_range"        : mf_shape_range,
            "mf_amp_range"          : mf_amp_range,
            "sweep_range"           : sweep_range,
            "pulse_time_threshold"  : pulse_time_threshold,
            "apulse_region"         : apulse_region,
            "resistance"            : resistance,
            "mf_shape_threshold"    : mf_shape_threshold,
            "mf_amp_threshold"      : mf_amp_threshold,
            "waveform_length"       : waveform_length,
            "baseline"              : baseline,
            "trigger_point"         : trigger_point,
            "integration"           : integration,
            "he_region"             : he_region
        }

        self.template_pmt_pulse = np.array([], dtype='float')

        self.histogram_names_list = ["pulse_charge_hist",
                                     "pulse_amplitude_hist",
                                     "apulse_charge_hist",
                                     "pulse_mf_shape_hist",
                                     "pulse_mf_amplitude_hist",
                                     "pulse_mf_shape_mf_amplitude_hist",
                                     "pulse_mf_shape_p_amplitude_hist",
                                     "pulse_mf_amplitude_p_amplitude",
                                     "pulse_peak_time_hist",
                                     "pulse_times_hist",
                                     "baseline_hist"
                                     ]

        self.histogram_dict = {}

    def set_up_histograms(self):
        pmt_pulse_charge_hist = ROOT.TH1F(self.get_pmt_id() + "_pulse_charge_spectrum",
                                          self.get_pmt_id() + "_pulse_charge_spectrum",
                                          self.get_nbins(),
                                          self.get_charge_range()[0],
                                          self.get_charge_range()[1])

        pmt_pulse_amplitude_hist = ROOT.TH1F(self.get_pmt_id() + "_pulse_amplitude_spectrum",
                                             self.get_pmt_id() + "_pulse_amplitude_spectrum",
                                             self.get_nbins(),
                                             self.get_amp_range()[0],
                                             self.get_amp_range()[1])

        pmt_apulse_charge_hist = ROOT.TH1F(self.get_pmt_id() + "_apulse_charge_spectrum",
                                           self.get_pmt_id() + "_apulse_charge_spectrum",
                                           self.get_nbins(),
                                           self.get_charge_range()[0],
                                           self.get_charge_range()[1])

        pmt_pulse_mf_shape_hist = ROOT.TH1F(self.get_pmt_id() + "_mf_shape",
                                            self.get_pmt_id() + "_mf_shape",
                                            self.get_nbins(),
                                            self.get_mf_shape_range()[0],
                                            self.get_mf_shape_range()[1])

        pmt_pulse_mf_amplitude_hist = ROOT.TH1F(self.get_pmt_id() + "_mf_amplitude",
                                                self.get_pmt_id() + "_mf_amplitude",
                                                self.get_nbins(),
                                                self.get_mf_amp_range()[0],
                                                self.get_mf_amp_range()[1])

        pmt_pulse_mf_shape_mf_amplitude_hist = ROOT.TH2F(self.get_pmt_id() + "_mf_shape_mf_amplitude",
                                                         self.get_pmt_id() + "_mf_shape_mf_amplitude",
                                                         self.get_nbins(),
                                                         self.get_mf_shape_range()[0],
                                                         self.get_mf_shape_range()[1],
                                                         self.get_nbins(),
                                                         self.get_mf_amp_range()[0],
                                                         self.get_mf_amp_range()[1])

        pmt_pulse_mf_shape_p_amplitude_hist = ROOT.TH2F(self.get_pmt_id() + "_mf_shape_p_amplitude",
                                                        self.get_pmt_id() + "_mf_shape_p_amplitude",
                                                        self.get_nbins(),
                                                        self.get_mf_shape_range()[0],
                                                        self.get_mf_shape_range()[1],
                                                        self.get_nbins(),
                                                        self.get_amp_range()[0],
                                                        self.get_amp_range()[1])

        pmt_pulse_mf_amplitude_p_amplitude_hist = ROOT.TH2F(self.get_pmt_id() + "_mf_amplitude_p_amplitude",
                                                            self.get_pmt_id() + "_mf_amplitude_p_amplitude",
                                                            self.get_nbins(),
                                                            self.get_mf_amp_range()[0],
                                                            self.get_mf_amp_range()[1],
                                                            self.get_nbins(),
                                                            self.get_amp_range()[0],
                                                            self.get_amp_range()[1])

        pmt_pulse_peak_time_hist = ROOT.TH1I(self.get_pmt_id() + "_pulse_peak_times",
                                             self.get_pmt_id() + "_pulse_peak_times",
                                             self.get_waveform_length(),
                                             0,
                                             self.get_waveform_length())

        pmt_pulse_times_hist = ROOT.TH1I(self.get_pmt_id() + "_pulse_times",
                                         self.get_pmt_id() + "_pulse_times",
                                         self.get_waveform_length(),
                                         0,
                                         self.get_waveform_length())

        pmt_baseline_hist = ROOT.TH1F(self.get_pmt_id() + "_baseline",
                                      self.get_pmt_id() + "_baseline",
                                      self.get_nbins(),
                                      self.get_baseline_setting() - 10,
                                      self.get_baseline_setting() + 10)

        self.histogram_dict = {
            "pulse_charge_hist": pmt_pulse_charge_hist,
            "pulse_amplitude_hist": pmt_pulse_amplitude_hist,
            "apulse_charge_hist": pmt_apulse_charge_hist,
            "pulse_mf_shape_hist": pmt_pulse_mf_shape_hist,
            "pulse_mf_amplitude_hist": pmt_pulse_mf_amplitude_hist,
            "pulse_mf_shape_mf_amplitude_hist": pmt_pulse_mf_shape_mf_amplitude_hist,
            "pulse_mf_shape_p_amplitude_hist": pmt_pulse_mf_shape_p_amplitude_hist,
            "pulse_mf_amplitude_p_amplitude_hist": pmt_pulse_mf_amplitude_p_amplitude_hist,
            "pulse_peak_time_hist": pmt_pulse_peak_time_hist,
            "pulse_time_hists": pmt_pulse_times_hist,
            "baseline_hist": pmt_baseline_hist
        }

    def get_template_bool(self):
        return self.template_bool

    def set_template_bool(self,  new_bool: bool):
        self.template_bool = new_bool

    @staticmethod
    def get_normalisation_factor(vector: list):
        norm = 0.0
        for i in range(len(vector)):
            norm += vector[i] * vector[i]
        return np.sqrt(norm)

    def get_setting_dict(self):
        return self.setting_dict

    def get_setting(self, description: str):
        return self.get_setting_dict()[description]

    def set_setting(self, description: str, value):
        if type(value) == list:
            assert len(value) <= 2
        else:
            pass
        self.setting_dict[description] = value

    def get_resistance(self):
        return self.get_setting("resistance")

    def set_resistance(self, new_resistance: float):
        self.set_setting("resistance", new_resistance)

    def get_trigger_point(self):
        return self.get_setting("trigger_point")

    def set_trigger_point(self, new_trigger_point: int):
        self.set_setting("trigger_point", new_trigger_point)

    def get_waveform_length(self):
        return self.get_setting("waveform_length")

    def set_waveform_length(self, new_waveform_length: int):
        self.set_setting("waveform_length", new_waveform_length)

    def get_apulse_region(self):
        return self.get_setting("apulse_region")

    def set_apulse_region(self, new_apulse_pos: int):
        self.set_setting("apulse_region", new_apulse_pos)

    def get_histogram_dict(self):
        return self.histogram_dict

    def get_histogram(self, description: str):
        return self.get_histogram_dict()[description]

    def get_sweep_range(self):
        return self.get_setting("sweep_range")

    def get_sweep_range_min(self):
        return self.get_setting("sweep_range")[0]

    def get_sweep_range_max(self):
        return self.get_setting("sweep_range")[1]

    def set_sweep_range(self, new_sweep_range: list):
        assert len(new_sweep_range) == 2
        assert new_sweep_range[1] > new_sweep_range[0]
        self.set_setting("sweep_range", new_sweep_range)

    def get_pulse_time_threshold(self):
        return self.get_setting("pulse_time_threshold")

    def set_pulse_time_threshold(self, new_threshold: int):
        self.set_setting("pulse_time_threshold", new_threshold)

    def get_nbins(self):
        return self.get_setting("nbins")

    def set_nbins(self, new_nbins: int):
        self.set_setting("nbins", new_nbins)

    def get_baseline_setting(self):
        return self.get_setting("baseline")

    def set_baseline_setting(self, new_baseline: int):
        self.set_setting("baseline", new_baseline)

    def get_charge_range(self):
        return self.get_setting("charge_range")

    def set_charge_range(self, new_charge_range: list):
        assert len(new_charge_range) == 2
        assert new_charge_range[1] > new_charge_range[0]
        self.set_setting("charge_range", new_charge_range)

    def get_amp_range(self):
        return self.get_setting("amp_range")

    def set_amp_range(self, new_amp_range: list):
        assert len(new_amp_range) == 2
        assert new_amp_range[1] > new_amp_range[0]
        self.set_setting("amp_range", new_amp_range)

    def get_mf_shape_range(self):
        return self.get_setting("mf_shape_range")

    def set_mf_shape_range(self, new_mf_shape_range: list):
        assert len(new_mf_shape_range) == 2
        assert new_mf_shape_range[1] > new_mf_shape_range[0]
        self.set_setting("mf_shape_range", new_mf_shape_range)

    def get_mf_amp_range(self):
        return self.get_setting("mf_amp_range")

    def set_mf_amp_range(self, new_mf_amp_range: list):
        assert len(new_mf_amp_range) == 2
        assert new_mf_amp_range[1] > new_mf_amp_range[0]
        self.set_setting("mf_amp_range", new_mf_amp_range)

    def get_histogram_names(self):
        return self.histogram_names_list

    def get_histogram_names_dict(self):
        return self.histogram_dict

    def get_pmt_id(self):
        return self.pmt_id

    def set_pmt_id(self, new_pmt_id: str):
        self.pmt_id = new_pmt_id

    def get_data_id(self):
        return self.data_id

    def get_charge_cut(self):
        return self.get_setting("charge_cut")

    def get_event_number(self):
        return self.number_of_events

    def set_number_of_events(self, new_event_number):
        self.number_of_events = new_event_number

    def create_pmt_pulse_template(self, root_file_name: str, template_histogram_name: str):
        template_root_file = ROOT.TFile(root_file_name, "READ")
        template_histogram = template_root_file.Get(template_histogram_name)
        template_list = []
        print('>>> Creating Template: ')
        print('>>> i \t V')
        for i_bin in range(int(template_histogram.GetEntries())):
            print(">>> ", i_bin, "\t", template_histogram.GetBinContent(i_bin))
            template_list.append(template_histogram.GetBinContent(i_bin))

        norm = self.get_normalisation_factor(template_list)

        self.set_template_pmt_pulse(np.array(template_list, dtype='float') / norm)
        # print(self.get_template_pmt_pulse())

        self.set_template_bool(True)

    def set_template_pmt_pulse(self, new_template_pmt_pulse: np.array):
        self.template_pmt_pulse = new_template_pmt_pulse

    def get_template_pmt_pulse(self):
        return self.template_pmt_pulse

    def get_pmt_pulse_charge_hist(self):
        return self.get_histogram("pulse_charge_hist")

    def fill_pmt_pulse_charge_hist(self, value: float):
        # print(value)
        self.get_histogram("pulse_charge_hist").Fill(value)

    def get_pmt_pulse_amplitude_hist(self):
        return self.get_histogram("pulse_amplitude_hist")

    def fill_pmt_pulse_amplitude_hist(self, value: float):
        self.get_histogram("pulse_amplitude_hist").Fill(value)

    def get_pmt_pulse_mf_shape_hist(self):
        return self.get_histogram("pulse_mf_shape_hist")

    def fill_pmt_pulse_mf_shape_hist(self, value: float):
        self.get_histogram("pulse_mf_shape_hist").Fill(value)

    def get_pmt_pulse_mf_amplitude_hist(self):
        return self.get_histogram("pulse_mf_amplitude_hist")

    def fill_pmt_pulse_mf_amplitude_hist(self, value: float):
        self.get_histogram("pulse_mf_amplitude_hist").Fill(value)

    def get_pmt_pulse_mf_shape_mf_amplitude(self):
        return self.get_histogram("pulse_mf_shape_mf_amplitude_hist")

    def fill_pmt_pulse_mf_shape_mf_amplitude(self, x_value: float, y_value: float):
        self.get_histogram("pulse_mf_shape_mf_amplitude_hist").Fill(x_value, y_value)

    def get_pmt_pulse_mf_shape_p_amplitude_hist(self):
        return self.get_histogram("pulse_mf_shape_p_amplitude_hist")

    def fill_pmt_pulse_mf_shape_p_amplitude_hist(self, x_value: float, y_value: float):
        self.get_histogram("pulse_mf_shape_p_amplitude_hist").Fill(x_value, y_value)

    def get_pmt_pulse_mf_amplitude_p_amplitude_hist(self):
        return self.get_histogram("pulse_mf_amplitude_p_amplitude_hist")

    def fill_pmt_pulse_mf_amplitude_p_amplitude_hist(self, x_value, y_value):
        self.get_histogram("pulse_mf_amplitude_p_amplitude_hist").Fill(x_value, y_value)

    def get_pmt_apulse_charge_hist(self):
        return self.get_histogram("apulse_charge_hist")

    def fill_pmt_apulse_charge_hist(self, value: float):
        self.get_histogram("apulse_charge_hist").Fill(value)

    def get_pmt_pulse_peak_time_hist(self):
        return self.get_histogram("pulse_peak_time_hist")

    def fill_pmt_pulse_peak_time_hist(self, value: int):
        self.get_histogram("pulse_peak_time_hist").Fill(value)

    def get_pmt_pulse_times_hist(self):
        return self.get_histogram("apulse_time_hists")

    def fill_pmt_pulse_times_hist(self, value: list):
        for i_value in range(len(value)):
            self.get_histogram("pulse_time_hists").Fill(value[i_value])

    def get_pmt_baseline_hist(self):
        return self.get_histogram("baseline_hist")

    def fill_pmt_baseline_hist(self, value: float):
        self.get_histogram("baseline_hist").Fill(value)

    def fill_pmt_hists(self, results: dict):

        # print(results)
        pulse_charge: float = results["pulse_charge"]
        pulse_amplitude: float = results["pulse_amplitude"]
        apulse_charge: float = results["apulse_charge"]
        mf_amplitude: float = results["pulse_mf_amp"]
        mf_shape: float = results["pulse_mf_shape"]
        pulse_peak_time: int = results["pulse_peak_time"]
        pulse_times: list = results["pulse_times"]
        baseline: float = results["baseline"]

        self.fill_pmt_pulse_charge_hist(pulse_charge)
        self.fill_pmt_pulse_amplitude_hist(pulse_amplitude)
        self.fill_pmt_apulse_charge_hist(apulse_charge)
        self.fill_pmt_pulse_mf_shape_hist(mf_shape)
        self.fill_pmt_pulse_mf_amplitude_hist(mf_amplitude)
        self.fill_pmt_pulse_mf_shape_mf_amplitude(mf_shape, mf_amplitude)
        self.fill_pmt_pulse_mf_shape_p_amplitude_hist(mf_shape, pulse_amplitude)
        self.fill_pmt_pulse_mf_amplitude_p_amplitude_hist(pulse_amplitude, mf_amplitude)
        self.fill_pmt_pulse_peak_time_hist(pulse_peak_time)
        self.fill_pmt_pulse_times_hist(pulse_times)
        self.fill_pmt_baseline_hist(baseline)

        self.set_number_of_events(self.get_event_number() + 1)

    def save_to_file(self, root_file_name: str):
        root_file = ROOT.TFile(root_file_name, "RECREATE")
        root_file.cd()
        for hist in self.get_histogram_dict().keys():
            self.get_histogram_dict()[hist].Write()
        root_file.Close()

    def save_histograms(self, directory: ROOT.TDirectory):
        directory.cd()
        for hist in self.get_histogram_dict().keys():
            self.get_histogram_dict()[hist].Write()

    @staticmethod
    def save_histogram(root_file: ROOT.TFile, hist, write_function: str):
        if write_function in ['RECREATE', "CREATE", "UPDATE"]:
            file = ROOT.TFile(root_file, write_function)
            file.cd()
            hist.Write()
        else:
            print("Invalid write function.")

    def set_sweep_bool(self, new_bool: bool):
        self.sweep_bool = new_bool

    def get_sweep_bool(self):
        return self.sweep_bool


class PMT_Waveform:

    def __init__(self, pmt_waveform_list: list, pmt_object: PMT_Object):

        # Store the PMT_Object into memory
        self.pmt_object = pmt_object

        self.pmt_trace_id = self.pmt_object.get_data_id() + "_ch" + self.pmt_object.get_pmt_id() + "_tr" + str(self.pmt_object.get_event_number())

        # If the input waveform list is of strings we need to check each element
        if type(pmt_waveform_list[0]) == str:
            temp_list_1 = []
            for i in range(len(pmt_waveform_list)):
                if pmt_waveform_list[i] == '' or pmt_waveform_list[i] == '\n':
                    pass
                else:
                    temp_list_1.append(pmt_waveform_list[i].strip())
        # If not a string we don't need to worry
        else:
            temp_list_1 = pmt_waveform_list

        # Store waveform
        self.pmt_waveform = np.array(temp_list_1, dtype='float')

        # Check whether the pulse lies where it is expected i.e above the 100 ns timestamp
        self.pmt_pulse_peak_amplitude = 0.0
        self.pmt_pulse_peak_position = np.argmin(self.pmt_waveform)

        self.pmt_pulse_start = 0
        self.pmt_pulse_end = 0
        self.pmt_baseline = 0.0
        self.pmt_pulse_charge = 0.0
        self.pmt_apulse_charge = 0.0
        self.pmt_pulse_mf_shape = 0.0
        self.pmt_pulse_mf_amp = 0.0
        self.pmt_pulse_trigger = False
        self.pmt_apulse_trigger = False
        self.pmt_pulse_times = []
        self.pmt_waveform_sweep_shape = []
        self.pmt_waveform_sweep_amp = []
        self.pmt_waveform_reduced = np.array([], dtype='float')
        self.done_sweep = False

        self.results_dict = {}
        if self.get_pmt_pulse_peak_position() < self.get_pmt_object().get_pulse_time_threshold():
            pass
        else:

            self.pmt_waveform_length = self.get_pmt_waveform().size
            self.set_pmt_baseline(np.average(self.pmt_waveform[0:self.get_pmt_object().get_setting("trigger_point")]))
            self.set_pmt_waveform_reduced()
            self.set_pmt_pulse_peak_amplitude(-1 * np.amin(self.get_pmt_waveform_reduced()))
            self.calculate_charge()

            if self.get_pmt_pulse_charge() < self.get_pmt_object().get_setting("charge_cut"):
                self.update_results_dict()
                return
            else:
                self.set_pmt_pulse_trigger(True)

                # Only sweep the waveform if there is a template
                if self.get_pmt_object().get_template_bool() and self.get_pmt_pulse_peak_position() < self.get_pmt_object().get_waveform_length() - self.get_pmt_object().get_template_pmt_pulse().size:

                    pmt_pulse = self.get_pmt_waveform()[self.get_pmt_pulse_peak_position() - np.argmin(self.get_pmt_object().get_template_pmt_pulse()): self.get_pmt_pulse_peak_position() - np.argmin(self.get_pmt_object().get_template_pmt_pulse()) + self.get_pmt_object().get_template_pmt_pulse().size] - self.get_pmt_baseline()

                    inner_product = np.dot(self.pmt_object.get_template_pmt_pulse(), pmt_pulse)
                    self.get_pmt_object().get_normalisation_factor(pmt_pulse)

                    self.set_pmt_pulse_mf_amp(inner_product)
                    self.set_pmt_pulse_mf_shape(inner_product/self.get_pmt_object().get_normalisation_factor(pmt_pulse))

                    if self.get_pmt_object().get_sweep_bool():
                        self.done_sweep = True
                        self.pmt_pulse_sweep()

        self.update_results_dict()
        # Store the pmt waveform length
        self.pmt_waveform_length = self.pmt_waveform.size
        # self.fill_pmt_hists()
        return

    def get_pulse_trigger(self):
        return self.pmt_pulse_trigger

    def get_pmt_object(self):
        return self.pmt_object

    def get_pmt_pulse_peak_amplitude(self):
        return self.pmt_pulse_peak_amplitude

    def set_pmt_pulse_peak_amplitude(self, new_amplitude: float):
        self.pmt_pulse_peak_amplitude = new_amplitude

    def get_pmt_apulse_charge(self):
        return self.pmt_apulse_charge

    def set_pmt_apulse_charge(self, new_apulse_charge: float):
        self.pmt_apulse_charge = new_apulse_charge

    def get_pmt_pulse_start(self):
        return self.pmt_pulse_start

    def get_results_dict(self):
        return self.results_dict

    def set_results_dict(self, description: str, value):
        self.get_results_dict()[description] = value

    def set_pmt_pulse_start(self, new_pmt_pulse_start: int):
        self.pmt_pulse_start = new_pmt_pulse_start

    def get_pmt_pulse_end(self):
        return self.pmt_pulse_end

    def set_pmt_pulse_end(self, new_pmt_pulse_end: int):
        self.pmt_pulse_end = new_pmt_pulse_end

    def get_pmt_pulse_trigger(self):
        return self.pmt_pulse_trigger

    def set_pmt_pulse_trigger(self, new_pmt_pulse_trigger: bool):
        self.pmt_pulse_trigger = new_pmt_pulse_trigger

    def get_pmt_waveform_reduced(self):
        return self.pmt_waveform_reduced

    def set_pmt_waveform_reduced(self):
        self.pmt_waveform_reduced = self.get_pmt_waveform() - self.get_pmt_baseline()

    def get_pmt_baseline(self):
        return self.pmt_baseline

    def set_pmt_baseline(self, new_pmt_baseline: float):
        self.pmt_baseline = new_pmt_baseline

    def get_pmt_pulse_charge(self):
        return self.pmt_pulse_charge

    def set_pmt_pulse_charge(self, new_pmt_pulse_charge: float):
        self.pmt_pulse_charge = new_pmt_pulse_charge

    def get_pmt_pulse_peak_position(self):
        return self.pmt_pulse_peak_position

    def get_pmt_waveform(self):
        return self.pmt_waveform

    def get_pmt_apulse_trigger(self):
        return self.pmt_apulse_trigger

    def set_pmt_apulse_trigger(self, new_bool: bool):
        self.pmt_apulse_trigger = new_bool

    def get_pmt_pulse(self):
        return self.pmt_waveform[self.get_pmt_pulse_start():self.get_pmt_pulse_end()] - self.pmt_baseline

    def get_pmt_waveform_length(self):
        return self.pmt_waveform_length

    def get_pmt_trace_id(self):
        return self.pmt_trace_id

    def fill_pmt_hists(self):
        #print(self.get_results_dict())
        self.pmt_object.fill_pmt_hists(self.get_results_dict())

    def update_results_dict(self):
        self.results_dict = {
            "pulse_charge": self.get_pmt_pulse_charge(),
            "pulse_amplitude": self.get_pmt_pulse_peak_amplitude(),
            "apulse_charge": self.get_pmt_apulse_charge(),
            "pulse_mf_shape": self.get_pmt_pulse_mf_shape(),
            "pulse_mf_amp": self.get_pmt_pulse_mf_amp(),
            "pulse_peak_time": self.get_pmt_pulse_peak_position(),
            "pulse_times": self.get_pmt_pulse_times(),
            "baseline": self.get_pmt_baseline()
        }

    def save_pmt_waveform_histogram(self, root_file: ROOT.TFile):
        name = self.pmt_object.get_pmt_id() + self.get_pmt_trace_id()
        pmt_waveform_hist = ROOT.TH1I(name, name, self.get_pmt_waveform_length(), 0, self.get_pmt_waveform_length())
        for timestamp in range(self.get_pmt_waveform_length()):
            pmt_waveform_hist.SetBinContent(timestamp, self.get_pmt_waveform()[timestamp])
        root_file.cd()
        pmt_waveform_hist.Write()
        del pmt_waveform_hist

    def check_cuts(self):
        # This function should be used in the init function and in the sweep function
        pass

    def pmt_pulse_sweep(self):
        sweep_start = self.get_pmt_object().get_sweep_range()[0]
        sweep_end = self.get_pmt_object().get_sweep_range()[1]
        sweep_window_length = self.pmt_object.get_template_pmt_pulse().size

        # plt.plot(self.get_pmt_waveform_reduced())
        # plt.xlabel('timestamp (ns)')
        # plt.title("Waveform {}".format(self.get_pmt_trace_id()))
        # plt.show(block=True)

        matched_filter_shape_list = []
        matched_filter_amplitude_list = []
        for i_sweep in range(sweep_start, sweep_end - sweep_window_length):
            pmt_waveform_section = self.get_pmt_waveform_reduced()[
                                   i_sweep : i_sweep + sweep_window_length]

            if pmt_waveform_section.size == self.pmt_object.get_template_pmt_pulse().size:
                pass
            else:
                print("Section size {} template size {}".format(pmt_waveform_section.size,
                                                                self.pmt_object.get_template_pmt_pulse().size))
                continue

            inner_product = np.dot(self.pmt_object.get_template_pmt_pulse(), pmt_waveform_section)

            matched_filter_amplitude_list.append(inner_product)
            matched_filter_shape_list.append(inner_product / np.sqrt(np.dot(pmt_waveform_section, pmt_waveform_section)))

            '''fig, ax1 = plt.subplots()
            color = 'tab:red'
            ax1.plot(pmt_waveform_section, color=color)
            ax1.set_xlabel('timestamp (ns)')
            ax1.set_ylabel('ADC /mV', color=color)
            ax1.tick_params(axis='y', labelcolor=color)
            color = 'tab:blue'
            ax2 = ax1.twinx()
            ax2.plot(self.pmt_object.get_template_pmt_pulse(), color=color)
            ax2.set_ylabel('Normalised ADC', color=color)  # we already handled the x-label with ax1
            ax2.tick_params(axis='y', labelcolor=color)
            fig.tight_layout()
            plt.text(0, 0, "Shape: {}".format(inner_product / np.sqrt(np.dot(pmt_waveform_section, pmt_waveform_section))))
            plt.show(block=False)
            plt.pause(0.01)
            plt.close()'''

        matched_filter_shape = np.array(matched_filter_shape_list)
        matched_filter_amplitude = np.array(matched_filter_amplitude_list)

        # TODO: check if both thresholds have been breached
        shape_peaks, _ = find_peaks(matched_filter_shape, height=self.get_pmt_object().get_setting("mf_shape_threshold"), distance=int(sweep_window_length / 2))
        amplitude_peaks, _ = find_peaks(matched_filter_amplitude, height=self.get_pmt_object().get_setting("mf_amp_threshold"), distance=int(sweep_window_length / 2))

        temp = []
        if len(shape_peaks) > 0:
            for index, value in enumerate(shape_peaks):
                if value in amplitude_peaks:
                    temp.append(value)

        if len(temp) > 0:
            self.set_pmt_apulse_trigger(True)
            self.set_pmt_pulse_times(temp)

        self.pmt_waveform_sweep_shape = matched_filter_shape
        self.pmt_waveform_sweep_amp = matched_filter_amplitude
        self.update_results_dict()

        '''fig, ax1 = plt.subplots()

        color = 'tab:red'
        ax1.set_xlabel('timestamp (ns)')
        ax1.set_ylabel('shape', color=color)
        ax1.plot(matched_filter_shape, color=color)
        ax1.plot(shape_peaks, matched_filter_shape[shape_peaks], "x", color='tab:green')
        ax1.tick_params(axis='y', labelcolor=color)
        ax1.plot(np.zeros_like(matched_filter_shape), "--", color="gray")

        ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

        color = 'tab:blue'
        ax2.set_ylabel('Waveform ADC /mV', color=color)  # we already handled the x-label with ax1
        ax2.plot(self.get_pmt_waveform_reduced(), color=color)
        ax2.tick_params(axis='y', labelcolor=color)

        ax2.set_ylabel('Amplitude', color=color)  # we already handled the x-label with ax1
        ax2.plot(matched_filter_amplitude, color=color)
        ax2.plot(amplitude_peaks, matched_filter_amplitude[amplitude_peaks], "x", color='tab:orange')
        ax2.tick_params(axis='y', labelcolor=color)

        fig.tight_layout()  # otherwise the right y-label is slightly clipped
        plt.show(block=True)'''

    def get_pmt_pulse_mf_shape(self):
        return self.pmt_pulse_mf_shape

    def set_pmt_pulse_mf_shape(self, new_shape: float):
        self.pmt_pulse_mf_shape = new_shape

    def get_pmt_pulse_mf_amp(self):
        return self.pmt_pulse_mf_amp

    def set_pmt_pulse_mf_amp(self, new_amp: float):
        self.pmt_pulse_mf_amp = new_amp

    def get_pmt_pulse_times(self):
        return self.pmt_pulse_times

    def set_pmt_pulse_times(self, new_pulse_times: list):
        self.pmt_pulse_times = new_pulse_times

    def calculate_charge(self):
        start = int(self.get_pmt_pulse_peak_position()) - 20
        end = int(self.get_pmt_pulse_peak_position()) + 30
        self.set_pmt_pulse_charge(-1 * (np.sum(self.get_pmt_waveform()[start:end] - self.get_pmt_baseline())) / self.get_pmt_object().get_resistance())
        '''for i in range(self.get_pmt_object().get_setting("trigger_point"), int(self.get_pmt_pulse_peak_position())):
            if self.get_pmt_waveform()[i] <= self.get_pmt_waveform()[self.get_pmt_pulse_peak_position()]/self.get_pmt_object().get_setting("integration")[0]:
                start = i
                break
        for i in range(int(self.get_pmt_pulse_peak_position()), self.get_pmt_waveform_length()):
            if self.get_pmt_waveform()[i] >= self.get_pmt_waveform()[self.get_pmt_pulse_peak_position()]/self.get_pmt_object().get_setting("integration")[1]:
                end = i
                break
        self.set_pmt_pulse_start(start)
        self.set_pmt_pulse_end(end)
        self.set_pmt_pulse_charge(-1 * (np.sum(self.get_pmt_waveform()[
                                               start:end] - self.get_pmt_baseline())) / self.get_pmt_object().get_resistance())'''


class PMT_Array:

    def __init__(self, topology: list, data_id: str):
        self.pmt_topology = topology
        self.pmt_object_array = []
        self.pmt_total_number = topology[0] * topology[1]

        for i_row in range(topology[0]):
            for i_col in range(topology[1]):
                pmt_number = i_col + i_row*topology[1]
                pmt_object = PMT_Object(str(pmt_number), data_id)
                self.append_pmt_object_array(pmt_object)
                del pmt_object

    def save_to_file(self, output_root_filename: str):
        output_root_file = ROOT.TFile(output_root_filename, "RECREATE")

        for i_pmt in range(self.get_pmt_total_number()):
            if self.get_pmt_object_number(i_pmt).get_event_number() > 0:
                output_root_file.cd()
                output_root_file.mkdir(self.get_pmt_object_number(i_pmt).get_pmt_id())
                directory = output_root_file.GetDirectory(self.get_pmt_object_number(i_pmt).get_pmt_id())
                self.get_pmt_object_number(i_pmt).save_histograms(directory)

    def get_pmt_topology(self):
        return self.pmt_topology

    def get_pmt_object_array(self):
        return self.pmt_object_array

    def get_pmt_total_number(self):
        return self.pmt_total_number

    def append_pmt_object_array(self, pmt_object: PMT_Object):
        self.pmt_object_array.append(pmt_object)

    def get_pmt_object_position(self, pmt_position: list):
        assert len(pmt_position) < 3
        if len(pmt_position) == 1:
            pmt_number = pmt_position[0]
        else:
            pmt_number = pmt_position[0] * self.get_pmt_topology()[1] + pmt_position[1]
        return self.pmt_object_array[pmt_number]

    def get_pmt_object_number(self, pmt_number: int):
        return self.get_pmt_object_array()[pmt_number]

    def set_pmt_templates(self, template_root_file_name: str, template_histogram_name_list: list):
        for i_pmt in range(self.get_pmt_total_number()):
            print(f">>> Setting template: {i_pmt}")
            self.get_pmt_object_number(i_pmt).create_pmt_pulse_template(template_root_file_name, template_histogram_name_list[i_pmt])

    def apply_setting(self, config_file_name: str):
        if config_file_name is not None:
            stuff_list = self.read_config_file(config_file_name)

            for i_pmt in range(self.get_pmt_total_number()):
                for i_setting in range(len(stuff_list)):
                    # print(self.get_pmt_object_number(i_pmt).get_setting_dict())
                    self.get_pmt_object_number(i_pmt).set_setting(stuff_list[i_setting][0], stuff_list[i_setting][1])
                    # print(self.get_pmt_object_number(i_pmt).get_setting_dict())
                self.get_pmt_object_number(i_pmt).set_up_histograms()
        else:
            for i_pmt in range(self.get_pmt_total_number()):
                self.get_pmt_object_number(i_pmt).set_up_histograms()

    def set_pmt_id(self, pmt_id: str, pmt_object_number: int):
        self.get_pmt_object_number(pmt_object_number).set_pmt_id(pmt_id)

    def set_sweep_bool(self, new_bool: bool):
        for i_pmt in range(self.get_pmt_total_number()):
            self.get_pmt_object_number(i_pmt).set_sweep_bool(new_bool)

    @staticmethod
    def read_config_file(config_file_name: str):
        try:
            config_file = open(config_file_name, 'r')
        except FileNotFoundError as fnf_error:
            print(fnf_error)
            raise Exception("Error opening config file")

        output_list = []
        for i_line, line in enumerate(config_file.readlines()):
            tokens = line.split(" ")
            if tokens[0] == '#':
                description = tokens[1].strip()
                value = tokens[3].split(",")
                value_ = []
                if len(value) == 1:
                    try:
                        value_ = int(value[0].strip())
                    except:
                        value_ = float(value[0].strip())
                else:
                    for i in range(len(value)):
                        try:
                            value_.append(int(value[i].strip()))
                        except:
                            value_.append(float(value[i].strip()))

                temp_list = [description, value_]
                output_list.append(temp_list)

        return output_list

    def fit_bismuth_function(self):
        # TODO: This function needs OOP

        for i_pmt in range(self.get_pmt_total_number()):
            if "GAO607" in self.get_pmt_object_number(i_pmt).get_pmt_id():
                if self.get_pmt_object_number(i_pmt).get_pmt_pulse_charge_hist().GetEntries() == 0:
                    pass
                else:
                    canvas_name = self.get_pmt_object_number(i_pmt).get_pmt_id()
                    canvas = ROOT.TCanvas(canvas_name, canvas_name)
                    fit = ROOT.TF1("fit",
                                   "[0]*"
                                   "(7.08*TMath::Gaus(x,[1],[2]) "
                                   " + 1.84*TMath::Gaus(x,[1]*(1 + 72.144/975.651),[2]*1.036) "
                                   " + 0.44*TMath::Gaus(x,[1]*(1 + 84.154/975.651),[2]*1.042)) "
                                   " + 0.464*(exp(0.254*x)/(1 + exp((x - 28.43)/2.14)))",
                                   33, 41)

                    fit.SetParNames("A", "mu", "sigma")

                    fit.SetParLimits(0, 0, 400)
                    fit.SetParLimits(1, 34, 37)
                    fit.SetParLimits(2, 0.8, 2)
                    fit.SetParameters(319, 36, 1.09)

                    self.get_pmt_object_number(i_pmt).get_pmt_pulse_charge_hist().Fit("fit", "", "", 33, 41)
                    self.get_pmt_object_number(i_pmt).get_pmt_pulse_charge_hist().SetXTitle("Charge /pC")
                    self.get_pmt_object_number(i_pmt).get_pmt_pulse_charge_hist().SetYTitle("Counts")
                    self.get_pmt_object_number(i_pmt).get_pmt_pulse_charge_hist().SetTitle("Bi Integrated Charge Spectrum")

                    resolution = (fit.GetParameter(2) / fit.GetParameter(1)) * 100.0
                    chi2 = fit.GetChisquare() / fit.GetNDF()
                    mu = fit.GetParameter(1)
                    mu_err = fit.GetParError(1)
                    sigma = fit.GetParameter(2)
                    sigma_err = fit.GetParError(2)

                    canvas.SetGrid()
                    canvas.Update()

                    canvas.Draw()
                    ROOT.gStyle.SetOptStat(11111111)
                    ROOT.gStyle.SetOptFit(1)
                    ROOT.gStyle.SetStatY(0.9)
                    ROOT.gStyle.SetStatX(0.9)
                    ROOT.gStyle.SetStatW(0.8)
                    ROOT.gStyle.SetStatH(0.1)
                    canvas.SaveAs(canvas_name+".pdf", "pdf")

            elif "GAO612" in self.get_pmt_object_number(i_pmt).get_pmt_id():
                if self.get_pmt_object_number(i_pmt).get_pmt_pulse_charge_hist().GetEntries() == 0:
                    pass
                else:
                    canvas_name = self.get_pmt_object_number(i_pmt).get_pmt_id()
                    canvas = ROOT.TCanvas(canvas_name, canvas_name)

                    fit = ROOT.TF1("fit",
                                   "[0]*"
                                   "(7.08*TMath::Gaus(x,[1],[2]) "
                                   " + 1.84*TMath::Gaus(x,[1]*(1 + 72.144/975.651),[2]*1.036) "
                                   " + 0.44*TMath::Gaus(x,[1]*(1 + 84.154/975.651),[2]*1.042)) "
                                   " + 0.515*(exp(0.2199*x)/(1 + exp((x - 31.68)/2.48)))",
                                   37, 44)

                    fit.SetParNames("A", "mu", "sigma")

                    fit.SetParLimits(0, 0, 400)
                    fit.SetParLimits(1, 37, 40)
                    fit.SetParLimits(2, 0.8, 2)

                    fit.SetParameters(319, 39, 1.09)

                    self.get_pmt_object_number(i_pmt).get_pmt_pulse_charge_hist().Fit("fit", "", "", 37, 44)
                    self.get_pmt_object_number(i_pmt).get_pmt_pulse_charge_hist().SetXTitle("Charge /pC")
                    self.get_pmt_object_number(i_pmt).get_pmt_pulse_charge_hist().SetYTitle("Counts")
                    self.get_pmt_object_number(i_pmt).get_pmt_pulse_charge_hist().SetTitle("Bi Integrated Charge Spectrum")
                    self.get_pmt_object_number(i_pmt).get_pmt_pulse_charge_hist().Draw()

                    resolution = (fit.GetParameter(2) / fit.GetParameter(1)) * 100.0
                    chi2 = fit.GetChisquare() / fit.GetNDF()
                    mu = fit.GetParameter(1)
                    mu_err = fit.GetParError(1)
                    sigma = fit.GetParameter(2)
                    sigma_err = fit.GetParError(2)

                    canvas.SetGrid()
                    canvas.Update()

                    canvas.Draw()
                    ROOT.gStyle.SetOptStat(11111111)
                    ROOT.gStyle.SetOptFit(1)
                    ROOT.gStyle.SetStatY(0.9)
                    ROOT.gStyle.SetStatX(0.9)
                    ROOT.gStyle.SetStatW(0.8)
                    ROOT.gStyle.SetStatH(0.1)
                    canvas.SaveAs(canvas_name+".pdf", "pdf")
