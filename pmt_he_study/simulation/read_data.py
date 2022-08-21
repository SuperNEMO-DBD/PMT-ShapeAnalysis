import sys
sys.path.insert(1, '../..')
from pmt_he_study.format_plot import *
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def main():
    output = {
        "injected_num": [],
        "injected_times": [],
        "injected_amp": [],
        "mf_times": [],
        "mf_shapes": [],
        "mf_amps": []
    }
    results = [[0 for i in range(100)], [0 for i in range(100)]]
    with open("/Users/williamquinn/Documents/PhD/PMT-ShapeAnalysis/pmt_he_study/simulation/data.txt", "r") as read_file:
        fl = read_file.readlines()
        for index, line in enumerate(fl):
            line_list = line.split(":")
            injected_times = [int(i.replace("]", "").replace("[", "")) for i in line_list[1].split(",")]
            injected_amps = [int(i.replace("]", "").replace("[", "")) for i in line_list[3].split(",")]
            mf_times = [int(i.replace("]", "").replace("[", "")) for i in line_list[5].split(",")]
            mf_shapes = [float(i.replace("]", "").replace("[", "")) for i in line_list[7].split(",")]
            mf_amps = [float(i.replace("]", "").replace("[", "")) for i in line_list[9].split(",")]

            selected = []
            rejected = []
            for ap_t in mf_times:
                try:
                    pos = injected_times.index(ap_t)
                    selected.append(pos)
                    continue
                except ValueError:
                    pass
                try:
                    pos = injected_times.index(ap_t + 1)
                    selected.append(pos)
                    continue
                except ValueError:
                    pass
                try:
                    pos = injected_times.index(ap_t - 1)
                    selected.append(pos)
                    continue
                except ValueError:
                    pass

            for i, ap_t in enumerate(injected_times):
                if ap_t in mf_times or ap_t-1 in mf_times or ap_t+1 in mf_times:
                    pass
                else:
                    rejected.append(i)

            injected_times = np.array(injected_times)
            injected_amps = np.array(injected_amps)

            output["injected_num"].append(len(injected_times))
            output["injected_times"].append(injected_times)
            output["injected_amp"].append(injected_amps)

            sel_times = injected_times[selected]
            sel_amps = injected_amps[selected]
            rej_times = injected_times[rejected]
            rej_amps = injected_amps[rejected]

            for amp in sel_amps:
                results[0][amp - 1] += 1
            for amp in rej_amps:
                results[1][amp - 1] += 1

            mf_times = np.array(mf_times)
            mf_shapes = np.array(mf_shapes)
            mf_amps = np.array(mf_amps)
            output["mf_times"].append(mf_times)
            output["mf_shapes"].append(mf_shapes)
            output["mf_amps"].append(mf_amps)

    df = pd.DataFrame(output)
    df.to_json("data.json")

    results[0] = np.array(results[0])
    results[1] = np.array(results[1])
    tot = results[0] + results[1]

    plt.figure(figsize=figsize)
    plt.bar([0.5 + i for i in range(1, 101)], results[0] / tot * 100, width=1, alpha=0.5, label='Success Rate')
    plt.bar([0.5 + i for i in range(1, 101)], results[1] / tot * 100, width=1, alpha=0.5, label='Failure Rate')
    plt.xlim(1, 101)
    plt.ylim(0, 100)
    plt.ylabel("Percentage")
    plt.xlabel("Injected Amplitude /mV")
    plt.tight_layout()
    plt.legend(loc='center right')
    plt.savefig("/Users/williamquinn/Desktop/success.pdf")


if __name__ == "__main__":
    main()
