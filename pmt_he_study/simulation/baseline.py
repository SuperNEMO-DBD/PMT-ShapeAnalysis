import numpy as np
import random
from tqdm import tqdm
from multiprocessing import Process, cpu_count
from time import time
import os


def generate_baseline_value(i, freq):
    a = 977
    pdf_lower = (freq / np.sum(freq))[0:6]
    pdf_higher = (freq / np.sum(freq))[6:]

    if i % 2 == 0:
        while 1:
            xi = random.randint(a, 980)
            yi = random.uniform(0, np.max(pdf_lower))

            if yi > pdf_lower[xi - 975]:
                continue
            else:
                break
        return float(xi)
    else:
        while 1:
            xi = random.randint(981, 983)
            yi = random.uniform(0, np.max(pdf_higher))

            if yi > pdf_higher[xi - 981]:
                continue
            else:
                break
        return float(xi)


def func(filename, freq, max_sims, num):
    with open(filename, 'w') as create_file:
        prev_time = time()
        for k in range(max_sims):
            waveform = []
            for i in range(7000):
                val = generate_baseline_value(i, freq)
                waveform.append(int(val))

            string = '<trace channel="0">'
            for i in range(7000):
                string += ' {}'.format(waveform[i])
            string += ' \n'
            create_file.write(string)

            interval = time() - prev_time
            prev_time = time()

            '''print(f"{num}: {interval}")
            print(f"{num}: {(interval * max_sims)/60}")
            break'''
            if (k+1) % 100 == 0:
                print('Number: {} {:.0}%'.format(num, (k+1)/max_sims * 100))


def main():
    baseline = []
    with open("/Users/williamquinn/Desktop/PMT_Project/res_data_ch0/Baseline.txt", "r") as file:
        fl = file.readlines()

        for i, line in enumerate(fl):
            line_list = line.split(" ")

            for j in range(len(line_list)):
                baseline.append(int(line_list[j].strip()))
    baseline = np.array(baseline)

    freq, bin_edges = np.histogram(baseline[:500], range=(975, 985), bins=10)

    max_apulse_num = 20
    max_sims = 1000 * max_apulse_num

    filename = "/Users/williamquinn/Documents/PhD/PMT-ShapeAnalysis/pmt_he_study/simulation/baselines.txt"
    create_file = open(filename, "w")
    for i in range(100):
        create_file.write("#blank \n")

    # n_cpu = cpu_count()
    n_cpu = 8

    sims = max_sims/n_cpu
    processes = [Process(target=func, args=(filename,
                                            freq, int(sims), t)) for t in range(n_cpu)]
    for p in processes:
        p.start()

    for p in processes:
        p.join()


if __name__ == "__main__":
    main()
