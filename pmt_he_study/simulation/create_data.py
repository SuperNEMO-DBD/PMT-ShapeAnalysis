import numpy as np
import ROOT


def sweep(vec: np.array, baseline: float, temp: np.array):
    temp_mf = {
        "apulse_amplitudes": [],
        "apulse_shapes": [],
        "apulse_times": [],
    }

    sweep_start = 800
    shape_cut = 0.8
    amp_cut = 1
    apulse_time_cut = 20

    apulse_amp_vec = []
    apulse_shape_vec = []
    apulse_time_vec = []

    current_apulse = 0
    previous_apulse = 0

    # shape_convolution = []
    # amp_convolution = []

    for i_sweep in range(sweep_start, len(vec) - temp.size):
        test = vec[i_sweep:i_sweep + temp.size] - baseline

        test_norm = np.sqrt(np.dot(test, test))
        amplitude_index = np.dot(test, temp)
        shape_index = amplitude_index / test_norm

        # print(shape_index)

        # shape_convolution.append(shape_index)
        # amp_convolution.append(amplitude_index)

        if shape_index > shape_cut and amplitude_index > amp_cut:
            distance_to_nearest_afterpulse = i_sweep - previous_apulse
            if distance_to_nearest_afterpulse > apulse_time_cut:
                apulse_amp_vec.append(amplitude_index)
                apulse_shape_vec.append(shape_index)
                apulse_time_vec.append(i_sweep)
                previous_apulse = i_sweep
                current_apulse = len(apulse_amp_vec) - 1
            else:
                if amplitude_index > apulse_amp_vec[current_apulse]:
                    apulse_amp_vec[current_apulse] = amplitude_index
                    apulse_shape_vec[current_apulse] = shape_index
                    apulse_time_vec[current_apulse] = i_sweep

    temp_mf['apulse_amplitudes'] = np.array(apulse_amp_vec)
    temp_mf['apulse_shapes'] = np.array(apulse_shape_vec)
    temp_mf['apulse_times'] = np.array(apulse_time_vec)
    return temp_mf


def get_baseline(waveform: np.array):
    temp = 0
    for i in range(500):
        temp += waveform[i]
    return temp / 500


def main():
    template_file = ROOT.TFile("/Users/williamquinn/Desktop/PMT_Project/pmt_short_templates.root")
    h_template = template_file.Get("Template_Ch0")
    x_temp = []
    y_temp = []

    for i in range(1, h_template.GetNbinsX() + 1):
        x_temp.append(i - 1)
        y_temp.append(h_template.GetBinContent(i))
    x_temp = np.array(x_temp, dtype='float')
    y_temp = np.array(y_temp, dtype='float')

    template = y_temp / np.sqrt(np.dot(y_temp, y_temp))

    apulse_times = []
    apulse_amplitudes = []
    data = []

    with open("/Users/williamquinn/Documents/PhD/PMT-ShapeAnalysis/pmt_he_study/simulation/apulses.txt", "r") as read_file:
        fl = read_file.readlines()
        for index, line in enumerate(fl):
            line_list = line.split(" ")
            if line_list[0] == "<ap_times:":
                times = line_list[1:101]
                for i in range(len(times)):
                    if i == 0:
                        times[i] = int(times[i].split("[")[-1].split(",")[0])
                    elif i == len(times) - 1:
                        times[i] = int(times[i].split("]")[0].split(",")[0])
                    else:
                        times[i] = int(times[i].split(",")[0])

                amps = line_list[102:]
                for i in range(len(amps)):
                    if i == 0:
                        amps[i] = int(amps[i].split("[")[-1].split(",")[0])
                    elif i == len(times) - 1:
                        amps[i] = int(amps[i].split("]")[0].split(",")[0])
                    else:
                        amps[i] = int(amps[i].split(",")[0])

                waveform = np.array(fl[index+1].split(" ")[2: -1], dtype=int)
                baseline = get_baseline(waveform)
                temp_mf = sweep(waveform, baseline, template)

                #apulse_times.append([times, temp_mf["apulse_times"]])
                #apulse_amplitudes.append([amps, temp_mf["apulse_amplitudes"]])

                data.append([times, amps, temp_mf])

    with open("/Users/williamquinn/Documents/PhD/PMT-ShapeAnalysis/pmt_he_study/simulation/data.txt", "w") as create_file:
        for i_w in range(len(data)):
            create_file.write("injected_times:{}:injected_amps:{}:mf_times:{}:mf_shapes:{}:mf_amps:{}:\n".format(data[i_w][0],
                                                                                                                 data[i_w][1],
                                                                                                                 list(data[i_w][2]["apulse_times"]),
                                                                                                                 list(data[i_w][2]["apulse_shapes"]),
                                                                                                                 list(data[i_w][2]["apulse_amplitudes"])))


if __name__ == "__main__":
    main()
