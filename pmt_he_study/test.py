import numpy as np
import scipy.fft as fft
import matplotlib.pyplot as plt
from scipy.signal import find_peaks


def parse_arguments():
    import argparse
    parser = argparse.ArgumentParser(description="Input file names")
    parser.add_argument('-i', required=True, type=str, help='Input data file')
    parser.add_argument('-t', required=True, type=str, help='template data file')
    args = parser.parse_args()
    return args


'''def read_file(filename: str):
    file_list = [[], []]
    file = open(filename)
    lines = file.readlines()
    for line in lines:
        line_list = line.strip(" ")
        print(line_list)
        file_list[0].append(int(line_list[0]))
        file_list[1].append(int(line_list[1]))
    file_list[0] = np.array(file_list[0])
    file_list[1] = np.array(file_list[1])
    return file_list'''


def normalise(temp: np.array):
    norm = np.dot(temp, temp)
    return temp/np.sqrt(norm)


def get_baseline(temp: np.array):
    return np.average(temp)


def convolution(data: np.array, template: np.array):
    conv = []
    template = normalise(template)
    size = template.size
    for i, i_data in enumerate(data[:-size]):
        test = data[i:i+size]
        test = normalise(test)

        temp_conv = np.dot(test, template)
        conv.append(temp_conv)
    conv = np.array(conv)
    # shape_peaks, _ = find_peaks(conv, height=0.9, distance=int(size/2))
    return conv


def fft_convolution(data: np.array, template: np.array):
    nfft = len(data)
    fft_data = fft.fft(data, nfft)
    fft_template = fft.fft(template, nfft)

    fft_conv = fft_data * np.conjugate(fft_template)
    conv = fft.ifft(fft_conv)
    # peaks, _ = find_peaks(conv, 5, distance=int(nfft/2))

    return conv


def plot_data(data: np.array):
    plt.plot(data[0], data[1])
    plt.show()
    plt.close()


def main():
    args = parse_arguments()
    input_file = args.i
    template_file = args.t

    data = np.loadtxt(input_file, unpack=True)
    template = np.loadtxt(template_file, unpack=True)

    conv = convolution(data[1], template[1])
    fft_peaks = fft_convolution(data[1], template[1])


if __name__ == '__main__':
    main()
