import sys

sys.path.insert(1, '../..')

import numpy as np
from scipy.optimize import curve_fit
import scipy.integrate as integrate
import csv
import ROOT
from array import array
from pmt_he_study.format_plot import *


def root_fit(E, y, yerr):
    x = pos(E)
    print(x, y, yerr)
    guess = array('d', [1E-15, pos(24.587), 1])
    # guess = array('d', [(y[-1] - y[-50])/((x[-1] - x[-50])*3600*24) /10000, y[0], 300*3600*24])
    names = ["A", "B", "C"]
    graph = ROOT.TGraphErrors(len(y), array('d', [i for i in x]), array('d', [i for i in y]),
                              array('d', [1 for i in range(len(x))]), array('d', [i for i in yerr]))

    fit = ROOT.TF1("func",
                   "[0]*(TMath::Log(x - [1]) + [2])/(x - [1])",
                   float(x[0]), float(x[-1]), len(guess))
    for i in range(len(guess)):
        fit.SetParameter(i, guess[i])
        fit.SetParName(i, names[i])
        # fit.SetParLimits(i, guess[i] - abs(guess[i])/2, guess[i] + abs(guess[i])/2)

    graph.Fit("func", "")
    pars = []
    errs = []
    chi = fit.GetChisquare()
    ndof = fit.GetNDF()

    for i in range(len(guess)):
        pars.append(fit.GetParameter(i))
        errs.append(fit.GetParError(i))
    return np.array(pars), np.array(errs), chi, ndof


def pos(E):
    V0 = 1400*0.42
    d = 10  # cm
    x = np.sqrt(E / V0) * d
    return x


def func(E, A, B, C):
    V0 = 1400*0.42
    d = 10  # cm
    x = pos(E)
    X = V0*(x/d)**2
    return A * (np.log(X - B) - C) / (X - B)


def func_x(x, A, B, C):
    V0 = 1400*0.42
    d = 10  # cm
    X = V0*(x/d)**2
    return A * (np.log(X - B) - C) / (X - B)


def func_e(E, A, B, C):
    return A * (np.log(E - B) - C) / (E - B)


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


def chi2(y_obs, y_exp, y_err, n_par):
    return np.sum(((y_obs - y_exp) / y_err) ** 2), len(y_obs) - n_par - 1


def main():
    energy = []
    xsection = []
    xsection_err = []

    with open('electron_xsections.csv') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                #print(f'Column names are {", ".join(row)}')
                line_count += 1
            else:
                energy.append(float(row[0]))
                xsection.append(float(row[1]) * 1e-17)
                xsection_err.append(float(row[2]) * 1e-17)
                #print(f'\t{row[0]} : {row[1]} ± {row[2]}')
                line_count += 1
        #print(f'Processed {line_count} lines.')

    energy = np.array(energy)
    position = np.array(pos(energy))
    xsection = np.array(xsection)
    xsection_err = np.array(xsection_err)

    double_energy = []
    double_xsection = []
    double_xsection_err = []

    with open('electron_double_xsections.csv') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                #print(f'Column names are {", ".join(row)}')
                line_count += 1
            else:
                double_energy.append(float(row[0]))
                double_xsection.append(float(row[1]) * 1e-19)
                double_xsection_err.append(float(row[2]) * 1e-19)
                #print(f'\t{row[0]} : {row[1]} ± {row[2]}')
                line_count += 1
        #print(f'Processed {line_count} lines.')

    double_energy = np.array(double_energy)
    double_position = np.array(pos(double_energy))
    double_xsection = np.array(double_xsection)
    double_xsection_err = np.array(double_xsection_err)

    guess = np.array([4E-15, -pos(25), -pos(3.93)])
    popt, pcov = curve_fit(f=func_e, xdata=energy, ydata=xsection, maxfev=100000)
    perr = np.sqrt(np.diag(pcov))

    x = np.array([x_i for index, x_i in enumerate(position) if 24.587 <= energy[index] <= 1400*0.42])
    e = np.array([e for e in energy if 24.587 <= e <= 1400*0.42])
    y = np.array([xsec for index, xsec in enumerate(xsection) if 24.587 <= energy[index] <= 1400*0.42])
    y_err = np.array([xsec_err for index, xsec_err in enumerate(xsection_err) if 24.587 <= energy[index] <= 1400*0.42])

    # popt, perr, chi, ndof = root_fit(e, y, y_err)
    model = func(e, *popt)
    chi, ndof = chi2(y, y_err, model, len(popt))

    fig1 = plt.figure(figsize=(5, 4), facecolor='w')
    frame1 = fig1.add_axes((.15, .34, .8, .6))
    frame1.set_xticklabels([])

    # plt.title(r'$e^-$ He Ionisation Cross Section vs Electron Energy')
    plt.errorbar(position, xsection, xsection_err, fmt='o', label="Single",
                 markersize=1, capsize=cap_size, linewidth=line_width, capthick=cap_thick, zorder=0)
    plt.errorbar(double_position, double_xsection, double_xsection_err, fmt='s', label="Double",
                 markersize=1, capsize=cap_size, linewidth=line_width, capthick=cap_thick, zorder=0)
    plt.plot(x, model, '-', label=r'Model $\chi^2$/N$_{DoF}:' + ' {:.2f}/{}$'.format(chi, ndof),
             linewidth=1, zorder=10, alpha=0.5)

    # plt.axvline(pos(1400*0.42), ls='--', color='black', linewidth=1.5, label='1st Dynode', zorder=20)
    # plt.xscale('log')
    plt.yscale('log')
    plt.ylabel(r'Cross Section /cm$^2$')
    plt.xlim(0, 10)
    plt.legend(loc="lower right")

    frame2 = fig1.add_axes((.15, .12, .8, .2))
    plt.errorbar(x, (y-model)/1E-18, yerr=y_err/1E-18, fmt="k.",
                 markersize=marker_size, capsize=cap_size, linewidth=line_width, capthick=cap_thick)
    plt.xlim(0, 10)
    # plt.xscale('log')
    plt.ylabel("model-data /1E-18")
    plt.xlabel('Displacement /cm')
    plt.axhline(0, ls='--', color='k')
    plt.tight_layout()

    print(popt, perr, chi, ndof)
    x0, x1 = pos(24.587), pos(1400*0.42)
    print(x0, x1)
    I = integrate.quad(func_x, x0, x1, args=(popt[0], popt[1], popt[2]))
    T = 292  # Kelvin
    k = 1.380649E-23  # Boltzman constant
    print("Integral:", I, "cm^3")
    print("err +   :", integrate.quad(func_x, x0, x1, args=(popt[0] + perr[0], popt[1] + perr[1], popt[2] + perr[2])))
    print("err -   :", integrate.quad(func_x, x0, x1, args=(popt[0] - perr[0], popt[1] - perr[1], popt[2] - perr[2])))
    print("Value   :", I[0]/1.380649E-23/292/(100*100*100))
    # print(func_integrand(energy, *popt, 24.587, 588))
    # I = func_integrand(energy, *popt, 24.587, 588)

    plt.savefig("/Users/williamquinn/Desktop/PMT_Project/xsec_energy.pdf")


if __name__ == "__main__":
    main()