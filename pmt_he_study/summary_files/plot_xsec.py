import sys

sys.path.insert(1, '../..')

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from functions.other_functions import process_date
import scipy.integrate as integrate
import csv
import ROOT
import time
from pmt_he_study.format_plot import *


def func(x, A, B, C):
    return A * (np.log(x - B) + C) / (x - B)


def func_integrand(x, A, B, C, x0, x1):
    return A * np.log((x1 - B) / (x0 - B)) * (0.5 * np.log((x1 - B) * (x0 - B)) + C)


def linear_voltage(r, V0, d):
    return V0 * r / d


def non_linear_voltage(r, V0, d):
    return V0 * (r / d) ** 2


def xsec(x, A, B):
    return (A * np.log(x) + B) / x


def poly(x, A, B, C, D):
    return A * x ** 3 + B * x ** 2 + C * x + D


def xsec_integrand(x, A, B):
    return (A * np.log(x) + B)


def chi2(y_obs, y_exp, y_err, Ndof):
    return np.sum(((y_obs - y_exp) / y_err) ** 2) / Ndof


def main():
    energy = []
    xsection = []
    xsection_err = []

    with open(
            '/Users/williamquinn/Documents/PhD/PMT-ShapeAnalysis/pmt_he_study/summary_files/electron_xsections.csv') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                print(f'Column names are {", ".join(row)}')
                line_count += 1
            else:
                energy.append(float(row[0]))
                xsection.append(float(row[1]) * 1e-17)
                xsection_err.append(float(row[2]) * 1e-17)
                print(f'\t{row[0]} : {row[1]} ± {row[2]}')
                line_count += 1
        print(f'Processed {line_count} lines.')

    energy = np.array(energy)
    xsection = np.array(xsection)
    xsection_err = np.array(xsection_err)

    double_energy = []
    double_xsection = []
    double_xsection_err = []

    with open(
            '/Users/williamquinn/Documents/PhD/PMT-ShapeAnalysis/pmt_he_study/summary_files/electron_double_xsections.csv') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                print(f'Column names are {", ".join(row)}')
                line_count += 1
            else:
                double_energy.append(float(row[0]))
                double_xsection.append(float(row[1]) * 1e-19)
                double_xsection_err.append(float(row[2]) * 1e-19)
                print(f'\t{row[0]} : {row[1]} ± {row[2]}')
                line_count += 1
        print(f'Processed {line_count} lines.')

    double_energy = np.array(double_energy)
    double_xsection = np.array(double_xsection)
    double_xsection_err = np.array(double_xsection_err)

    popt, pcov = curve_fit(f=func, xdata=energy, ydata=xsection, maxfev=10000)
    plt.figure(figsize=figsize, facecolor='white')
    plt.errorbar(energy, xsection, xsection_err, fmt='o', label="Single", markersize=marker_size, zorder=0)
    plt.errorbar(double_energy, double_xsection, double_xsection_err, fmt='s', label="Double",
                 markersize=marker_size, zorder=0)
    plt.plot(energy, func(energy, *popt), '-', label='Model', linewidth=1.5, zorder=10)
    plt.axvline(1400*0.42, ls='--', color='black', linewidth=1.5, label='1st Dynode', zorder=20)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Electron Energy /eV')
    plt.ylabel(r'Cross Section /$cm^2$')
    plt.title(r'$e^-$ He Ionisation Cross Section vs Electron Energy')
    plt.legend(loc='best')
    plt.tight_layout()

    print(popt)
    print(integrate.quad(func, 24.587, 1400*0.42, args=(popt[0], popt[1], popt[2])))
    print(func_integrand(energy, *popt, 24.587, 588))
    # I = func_integrand(energy, *popt, 24.587, 588)

    plt.savefig("/Users/williamquinn/Desktop/PMT_Project/xsec_energy.pdf")


if __name__ == "__main__":
    main()