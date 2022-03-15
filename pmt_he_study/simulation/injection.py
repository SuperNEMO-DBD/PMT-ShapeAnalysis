import numpy as np
import random
from multiprocessing import Process, cpu_count
import ROOT


def main():
    baselines = []
    with open("/Users/williamquinn/Documents/PhD/PMT-ShapeAnalysis/pmt_he_study/simulation/baselines.txt", "r") as read_file:
        fl = read_file.readlines()
        for index, line in enumerate(fl):
            line_list = line.split(" ")
            if len(line_list) < 7000:
                continue
            baselines.append([int(i) for i in line_list[2:-1]])

    template_file = ROOT.TFile("/Users/williamquinn/Desktop/PMT_Project/pmt_short_templates.root")
    h_template = template_file.Get("Template_Ch0")
    x_temp = []
    y_temp = []

    for i in range(1, h_template.GetNbinsX() + 1):
        x_temp.append(i - 1)
        y_temp.append(h_template.GetBinContent(i))
    x_temp = np.array(x_temp, dtype='float')
    y_temp = np.array(y_temp, dtype='float')

    template = -1 * y_temp / np.amin(y_temp)
    apulse_times = []
    apulse_amplitudes = []

    template_size = template.size

    for i in range(len(baselines)):
        ts = []
        amps = []
        for apn in range(100):
            time_ok = False
            ti = 0
            while not time_ok:
                ti = random.randint(800, 7000 - template_size)
                time_ok = True
                for j in range(len(ts)):
                    if ts[j] < ti < ts[j] + template_size:
                        time_ok = False
                    else:
                        pass

            ai = random.randrange(1, 101)
            baselines[i][ti:ti + template_size] += np.array(template * ai, dtype='int')

            ts.append(ti)
            amps.append(ai)

        apulse_times.append(ts)
        apulse_amplitudes.append(amps)

    with open("/Users/williamquinn/Documents/PhD/PMT-ShapeAnalysis/pmt_he_study/simulation/apulses.txt", 'w') as create_file:
        for i in range(100):
            create_file.write("#blank \n")
        for i in range(len(baselines)):
            string = f'<ap_times: {apulse_times[i]} ap_amps: {apulse_amplitudes[i]}>\n'
            create_file.write(string)
            string = '<trace channel="0">'
            for j in range(len(baselines[i])):
                string += ' {}'.format(baselines[i][j])
            string += ' \n'
            create_file.write(string)


if __name__ == "__main__":
    main()
