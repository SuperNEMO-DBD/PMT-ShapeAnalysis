import sys

import numpy as np

sys.path.insert(1, '../..')

import ROOT
from matplotlib.colors import LogNorm
from scipy.optimize import curve_fit
from pmt_he_study.format_plot import *
from functions.other_functions import *
from src.PMT_Classes import *
from mpl_toolkits.axes_grid1 import make_axes_locatable


def lorentzian(x, amp1, cen1, wid1):
    return amp1 * wid1 ** 2 / ((x - cen1) ** 2 + wid1 ** 2)


def main():
    topology = [2, 1]
    pmt_array = PMT_Array(topology, "summary")
    pmt_array.set_pmt_id("GAO607", 0)
    pmt_array.set_pmt_id("GAO612", 1)

    filename = "/Users/williamquinn/Desktop/data/set_5/200706_A1400_B1400_t0012.root"
    file = ROOT.TFile(filename, "READ")
    tree = file.T

    pulse_charges = []
    ap_region_charges = []
    he_ap_region_charges = []
    ratio = []
    he_ratio = []

    for event in tree:
        if event.OM_ID == 1:
            continue
        p_charge = event.pulse_charge
        ap_region_charge = event.ap_region_charge
        he_ap_region_charge = event.he_ap_region_charge

        pulse_charges.append(p_charge)
        ap_region_charges.append(ap_region_charge)
        he_ap_region_charges.append(he_ap_region_charge)

        ratio.append(ap_region_charge/p_charge)
        he_ratio.append(he_ap_region_charge/p_charge)

    plt.figure(figsize=figsize)
    freq, bin_edges = np.histogram(pulse_charges, range=(0, 400), bins=240)
    width = (bin_edges[2] - bin_edges[1])
    bin_centres = bin_edges[:-1] + width/2
    plt.bar(bin_centres, freq/np.sum(freq) * 100, width=width)
    plt.xlabel('Charge /pC')
    plt.ylabel('Relative Counts')
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/charge.pdf")

    plt.figure(figsize=figsize)
    freq_0, bin_edges = np.histogram(ap_region_charges, range=(0, 400), bins=240)
    width = (bin_edges[2] - bin_edges[1])
    bin_centres_0 = bin_edges[:-1] + width / 2
    plt.bar(bin_centres_0, freq_0 / np.sum(freq_0) * 100, width=width)
    plt.xlabel('Charge /pC')
    plt.ylabel('Relative Counts')
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/ap_charge.pdf")

    plt.figure(figsize=figsize)
    freq_1, bin_edges = np.histogram(he_ap_region_charges, range=(0, 400), bins=240)
    width = (bin_edges[2] - bin_edges[1])
    bin_centres_1 = bin_edges[:-1] + width / 2
    plt.bar(bin_centres_1, freq_1 / np.sum(freq_1) * 100, width=width)
    plt.xlabel('Charge /pC')
    plt.ylabel('Relative Counts')
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/he_ap_charge.pdf")

    plt.figure(figsize=(3, 2))
    freq_2, bin_edges = np.histogram(ratio, range=(-1, 1), bins=100)
    width = (bin_edges[2] - bin_edges[1])
    bin_centres_2 = bin_edges[:-1] + width / 2
    plt.errorbar(bin_centres_2, freq_2, yerr=np.sqrt(freq_2), fmt='k.', label='data', markersize=1,
                 capsize=0.5, linewidth=0.5, capthick=0.5)
    plt.bar(bin_centres_2, freq_2, width=width, alpha=0.5)
    plt.xlabel('Ratio')
    plt.ylabel('Counts')
    plt.xlim(-1, 1)
    popt, pcov = curve_fit(lorentzian, bin_centres_2, freq_2,
                           p0=[10000, 0, 5],
                           bounds=[[100, -0.25, 0], [100000, 0.25, 10]], maxfev=5000)
    plt.plot(np.linspace(-1, 1, 1000), lorentzian(np.linspace(-1, 1, 1000), *popt), "C2--", label='f(x) = Lorentz',
             linewidth=1)

    handles, labels = plt.gca().get_legend_handles_labels()
    '''patch = patches.Patch(color='white', label=r'$A =$ {:.3e} ± {:.0e}'.format(popt[0], np.sqrt(pcov[0, 0])))
    patch_1 = patches.Patch(color='white', label=r'$\mu =$ {:.4f} ± {:.4f}'.format(popt[1], np.sqrt(pcov[1, 1])))
    patch_2 = patches.Patch(color='white', label=r'$\sigma =$ {:.4f} ± {:.4f}'.format(popt[2], np.sqrt(pcov[2, 2])))
    patch_3 = patches.Patch(color='white', label=r'$\chi^2_R =$ {:.2f}'.format(chi2(freq_2, np.sqrt(freq_2), lorentzian(bin_centres_2, *popt), len(popt))))
    handles.extend([patch, patch_1, patch_2, patch_3])'''
    plt.legend(handles=handles, loc='best')
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/ratio.pdf")

    plt.figure(figsize=figsize)
    freq_3, bin_edges = np.histogram(he_ratio, range=(-1, 1), bins=240)
    width = (bin_edges[2] - bin_edges[1])
    bin_centres_3 = bin_edges[:-1] + width / 2
    plt.bar(bin_centres_3, freq_3 / np.sum(freq_3) * 100, width=width)
    plt.xlabel('Ratio')
    plt.ylabel('Relative Counts')
    plt.tight_layout()
    plt.axvline(np.average(he_ratio), ls='--', color='C1')
    plt.axvline(np.average(ap_region_charges) / np.average(pulse_charges), ls='--', color='C2')
    plt.savefig("/Users/williamquinn/Desktop/he_ratio.pdf")

    # print(np.std(ratio))

    fig = plt.figure(figsize=(3, 2))
    x = np.linspace(0, 400, 2)
    x_bins = np.linspace(0, 400, 50)
    y_bins = np.linspace(-50, 75, 50)
    plt.hist2d(pulse_charges, ap_region_charges, bins=[x_bins, y_bins], cmap=plt.get_cmap('Purples'), norm=LogNorm())
    plt.plot(x, x*np.average(ratio), '--', label='Ratio: {:.2f}'.format(np.average(ratio)),)
    plt.grid()
    plt.xlabel('Pulse Charge /pC')
    plt.ylabel('AP Region Charge /pC')
    plt.xlim(0, 400)
    plt.ylim(-50, 75)
    # plt.legend(loc='best')
    plt.colorbar()

    '''axs[1].bar(bin_centres_2, freq_2 / np.sum(freq_2) * 100, width=width)
    axs[1].set_xlabel('Ratio')
    axs[1].set_ylabel('Relative Counts')'''

    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/pc_vs_aprc.pdf")


if __name__ == "__main__":
    main()

