import sys

import matplotlib.pyplot as plt
import pandas as pd
from pandas.plotting import lag_plot, autocorrelation_plot
from sklearn.metrics import mean_squared_error
from statsmodels.tsa.ar_model import AutoReg
from statsmodels.tsa.ar_model import AR
from statsmodels.tsa.seasonal import seasonal_decompose
from sklearn.metrics import r2_score

sys.path.insert(1, '../..')

from pmt_he_study.models import *
from tqdm import tqdm

# import custom made classes
from functions.other_functions import do_bi_1000, chi2, process_date, linear, gaussian
from src.PMT_Array import PMT_Array


# bismuth_func_string_0 = "[0]*(7.08*TMath::Gaus(x,[1],[2]) + 1.84*TMath::Gaus(x,[1]*(1 + 72.144/975.651),[2]*1.036) + 0.44*TMath::Gaus(x,[1]*(1 + 84.154/975.651),[2]*1.042)) + 0.464*(exp(0.254*x)/(1 + exp((x - 28.43)/2.14)))"
# bismuth_func_string_1 = "[0]*(7.08*TMath::Gaus(x,[1],[2]) + 1.84*TMath::Gaus(x,[1]*(1 + 72.144/975.651),[2]*1.036) + 0.44*TMath::Gaus(x,[1]*(1 + 84.154/975.651),[2]*1.042)) + 0.515*(exp(0.2199*x)/(1 + exp((x - 31.68)/2.48)))"


def tot_bi(x: np.array, mu: float, sig: float, A: float):
    a = np.array(gaussian_noh(x, mu, sig, A * 7.08))
    b = np.array(gaussian_noh(x, mu * (1 + 72.144 / 975.651), sig * 1.036, A * 1.84))
    c = np.array(gaussian_noh(x, mu * (1 + 84.154 / 975.651), sig * 1.042, A * 0.44))
    return a + b + c


def plot_fit(charges: np.array, fit_pars, date, name):
    max_charge = 60
    min_charge = 0
    n_bins = 240

    popt = np.array([fit_pars[0][0], fit_pars[1][0], fit_pars[2][0]])
    mu, sig, A = popt[0], popt[1], popt[2]
    perr = [fit_pars[0][1], fit_pars[1][1], fit_pars[2][1]]
    mu_err, sig_err, A_err = perr[0], perr[1], perr[2]
    chi = fit_pars[-1][0]
    ndof = fit_pars[-1][1]

    lower = mu - 1
    higher = mu + 6

    freq, bin_centres, bin_width = charges[0], charges[1], charges[2]

    low_bin = int(lower / bin_width)
    high_bin = int(higher / bin_width)

    x = bin_centres[low_bin:high_bin]
    y = freq[low_bin:high_bin]

    fig1 = plt.figure(figsize=(5, 3), facecolor='white')
    frame1 = fig1.add_axes((.125, .3, .55, .6))

    plt.bar(bin_centres, freq, width=bin_width, color='b', alpha=0.25)
    plt.errorbar(x, y, yerr=np.sqrt(y), fmt='k.', label='Data', markersize=3, capsize=1, linewidth=1, capthick=1)
    plt.plot(x, tot_bi(x, *popt), 'r-', label='Model')
    plt.plot(x, gaussian_noh(x, mu, sig, A * 7.08), '--', label='CE K 976 keV')
    plt.plot(x, gaussian_noh(x, mu * (1 + 72.144 / 975.651), sig * 1.036, A * 1.84), '--', label='CE L 1048 keV')
    plt.plot(x, gaussian_noh(x, mu * (1 + 84.154 / 975.651), sig * 1.042, A * 0.44), '--', label='CE M 1060 keV')

    mcs = [r'$\mu$', r'$\sigma$', r'$A$']
    n_sf = get_n_sfs(popt, perr)
    f_strs = ['' for i in range(len(perr))]
    for j in range(len(n_sf)):
        f_strs[j] += mcs[j] + "= {:." + str(n_sf[j]) + "f}±{:." + str(n_sf[j]) + "f} "

    frame1.set_xticklabels([])
    handles, labels = plt.gca().get_legend_handles_labels()
    patch = patches.Patch(color='white', label=f_strs[0].format(mu, mu_err))
    patch_1 = patches.Patch(color='white', label=f_strs[1].format(sig, sig_err))
    patch_2 = patches.Patch(color='white', label=f_strs[2].format(A, A_err))
    patch_3 = patches.Patch(color='white', label=r'$\chi/N_{DoF} =$' + ' {:.2f}/{}'.format(chi, ndof))
    handles.extend([patch, patch_1, patch_2, patch_3])
    plt.legend(handles=handles, bbox_to_anchor=(1.0, 1))

    plt.ylabel('counts')
    plt.xlim(bin_centres[low_bin - 5], bin_centres[high_bin + 5])
    plt.title(date + ' Resolution @1MeV: {:.2f} %'.format(sig / mu * 100))

    # Residual plot
    frame2 = fig1.add_axes((.125, .15, .55, .125))
    temp = (tot_bi(x, *popt) - y) / tot_bi(x, *popt)
    temp_err = np.sqrt(y) / tot_bi(x, *popt)
    plt.errorbar(x, temp, yerr=temp_err, fmt='k.', markersize=3, capsize=1, linewidth=1, capthick=1)
    plt.xlim(bin_centres[low_bin - 5], bin_centres[high_bin + 5])
    plt.xlabel('charge /pC')
    plt.ylim(-1, 1)
    plt.ylabel('model - data /model', fontsize=6)
    plt.axhline(0, ls='--', color='black')
    plt.savefig("/Users/williamquinn/Desktop/PMT_Project/res_plots/" + date + "_1kV_fit_bi_spec_" + name + ".pdf")
    plt.close(fig1)


def read_file(date: str, voltage: int, root_file_name: str, pmt_array: PMT_Array):
    file = ROOT.TFile(root_file_name, "READ")
    file.cd()

    fit_parameter = [None for i in range(pmt_array.get_pmt_total_number())]

    for i_om in range(pmt_array.get_pmt_total_number()):
        charge_hist = file.Get(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() +
                               "_charge_spectrum_" + str(voltage) + "V")
        amp_hist = file.Get(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() +
                            "_amplitude_spectrum_" + str(voltage) + "V")
        baseline_hist = file.Get(date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() +
                                 "_baseline_distribution_" + str(voltage) + "V")

        try:
            charge_hist.GetEntries()
            amp_hist.GetEntries()
            baseline_hist.GetEntries()
        except AttributeError:
            # print(root_file_name, ": No entries :", date + "_" + pmt_array.get_pmt_object_number(i_om).get_pmt_id() + "_charge_spectrum_" + str(voltage) + "V")
            continue

        mu_guess = charge_hist.GetMaximumBin() * charge_hist.GetBinWidth(0)
        fit_vals = do_bi_1000(charge_hist, mu_guess)

        freq, bin_centres = [], []
        for i_bin in range(0, charge_hist.GetNbinsX()):
            freq.append(charge_hist.GetBinContent(i_bin + 1))
            bin_centres.append(i_bin * charge_hist.GetBinWidth(i_bin + 1) + charge_hist.GetBinWidth(i_bin + 1))
        freq, bin_centres, bin_width = np.array(freq), np.array(bin_centres), charge_hist.GetBinWidth(3)
        plot_fit([freq, bin_centres, bin_width], fit_vals, date, str(i_om))

        if fit_vals[-1][0] == -1:
            pass
        else:
            pars = {
                # TODO: add the correlation matrix to this
                "mu": fit_vals[0][0],
                "mu_err": fit_vals[0][1],
                "sig": fit_vals[1][0],
                "sig_err": fit_vals[1][1],
                "base_mu": baseline_hist.GetMean(),
                "base_sig": baseline_hist.GetStdDev(),
                "chi2": fit_vals[-1][0],
                "ndof": fit_vals[-1][1]
            }
            fit_parameter[i_om] = pars

    file.Close()
    return fit_parameter


def fit_straight_line(x, y, dy, guess):
    scale = abs(np.max(x) - np.min(x)) / abs(np.max(y) - np.min(y))

    X = x - np.average(x)
    Y = y * scale
    dY = dy * scale

    popt, pcov = curve_fit(f=linear, xdata=X, ydata=Y, sigma=dY, absolute_sigma=True, p0=guess)
    perr = np.sqrt(np.diag(pcov))
    d_inv = np.diag(1 / np.sqrt(np.diag(pcov)))
    pcor = d_inv @ pcov @ d_inv
    popt_ = np.array([popt[0] / scale, (popt[1] - popt[0] * np.average(x)) / scale])
    perr_ = np.array([perr[0] / scale, np.sqrt(perr[1] ** 2 + (np.average(x) * perr[0]) ** 2) / scale])

    fit = {"popt": popt_, "perr": perr_, "pcor": pcor}

    return fit


def plot_res():
    data = pd.read_csv("/Users/williamquinn/Desktop/PMT_Project/res_vs_time.csv", header=0,
                       dtype={0: int, 1: int, 2: float, 3: float, 4: float, 5: int, 6: float, 7: float},
                       engine='python')
    data["time"] = process_date(data["date"])
    data = data[(data['date'] == 191023) | (data['date'] == 200212) | (data['date'] == 200429)
                | (data['date'] == 200717) | (data['date'] == 201015) | (data['date'] == 201216)]

    plt.figure(figsize=figsize, facecolor='white')

    sel_dates = [191023, 200212, 200429, 200717, 201015, 201216]
    chi_cut = 10
    data_ch0 = data[data["pmt"] == 0]
    x, y, dy = data_ch0[data_ch0['chi_2'] / data_ch0['ndof'] < chi_cut]['time'].values, \
               data_ch0[data_ch0['chi_2'] / data_ch0['ndof'] < chi_cut]['res'].values * 100, \
               data_ch0[data_ch0['chi_2'] / data_ch0['ndof'] < chi_cut]['res_err'].values * 100

    exposed_fit = fit_straight_line(x, y, dy, guess=[0, 0])
    plt.errorbar(x, y, yerr=dy, fmt='k.', label='Exposed', markersize=3, capsize=1, linewidth=1, capthick=1)
    plt.plot([x[0], x[-1]], linear([x[0], x[-1]], *exposed_fit["popt"]), 'r-')

    chi_cut = 10
    data_ch1 = data[data["pmt"] == 1]
    x, y, dy = data_ch1[data_ch1['chi_2'] / data_ch1['ndof'] < chi_cut]['time'].values, \
               data_ch1[data_ch1['chi_2'] / data_ch1['ndof'] < chi_cut]['res'].values * 100, \
               data_ch1[data_ch1['chi_2'] / data_ch1['ndof'] < chi_cut]['res_err'].values * 100
    control_fit = fit_straight_line(x, y, dy, guess=[0, 0])
    plt.errorbar(x, y, yerr=dy, fmt='.', label='Control', markersize=3, capsize=1, linewidth=1, capthick=1,
                 color='grey', ecolor='grey')
    plt.plot([x[0], x[-1]], linear([x[0], x[-1]], *control_fit["popt"]), 'r-')

    plt.ylim(1.5, 5)

    plt.xlabel('Days from 1% He Onset')
    plt.ylabel('Resolution at 1MeV /%')
    plt.fill_between([-100, 0], [5, 5], alpha=0.1,
                     facecolor='green', label='Atmospheric He')
    plt.fill_between([0, 98], [5, 5], alpha=0.1,
                     facecolor='blue', label='1% He')
    plt.fill_between([98, 500], [5, 5], alpha=0.1,
                     facecolor='red', label='10% He')
    plt.xlim(-30, 430)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/PMT_Project/res.pdf")
    plt.close()


def store_res(pmt_array: PMT_Array):
    out_file = open(f"/Users/williamquinn/Desktop/PMT_Project/res_vs_time.csv", "w")
    out_file.write("pmt,date,res,res_err,chi_2,ndof,base_mean,base_sig\n")

    filenames_txt = "/Users/williamquinn/Desktop/PMT_Project/filenames.txt"

    try:
        print(">>> Reading data from file: {}".format(filenames_txt))
        date_file = open(filenames_txt, 'r')
    except FileNotFoundError:
        raise Exception("Error opening data file {}".format(filenames_txt))

    filenames = np.loadtxt(filenames_txt, delimiter=',', dtype={'names': ['filename'], 'formats': ['S100']},
                           unpack=True)

    for i_file in tqdm(range(filenames.size)):
        file = filenames[i_file][0].decode("utf-8")

        date = file.split("_")[0]
        voltage = int(file.split("_")[1].split("A")[1])

        if voltage != 1000:
            continue

        fit_pars = read_file(date, voltage, "/Users/williamquinn/Desktop/PMT_Project/data/1000V/" + file, pmt_array)

        for i_om in range(pmt_array.get_pmt_total_number()):

            if fit_pars[i_om] is None:
                continue

            mu = fit_pars[i_om]["mu"]
            mu_err = fit_pars[i_om]["mu_err"]
            sig = fit_pars[i_om]["sig"]
            sig_err = fit_pars[i_om]["sig_err"]
            chi_2 = fit_pars[i_om]["chi2"]
            ndof = fit_pars[i_om]["ndof"]
            baseline_mean = fit_pars[i_om]["base_mu"]
            baseline_sig = fit_pars[i_om]["base_sig"]

            res = sig / mu
            res_err = res * ((sig_err / sig) ** 2 + (mu_err / mu) ** 2)

            out_file.write(f'{i_om},{date},{res},{res_err},{chi_2},{ndof},{baseline_mean},{baseline_sig}\n')

    out_file.close()

    print("<<<< FINISHED >>>")


def plot_base_drift():
    data = pd.read_csv("/Users/williamquinn/Desktop/PMT_Project/res_vs_time.csv", header=0,
                        dtype={0: int, 1: int, 2: float, 3: float, 4: float, 5: int, 6: float, 7: float},
                        engine='python')
    data["time"] = process_date(data["date"])
    data = data[data["base_mean"] > 900]

    fig1 = plt.figure(figsize=figsize, facecolor='white')
    #frame1 = fig1.add_axes((.125, .3, .55, .6))

    data_ch0 = data[data["pmt"] == 0]
    exposed_fit = fit_straight_line(data_ch0['time'].values, data_ch0['base_mean'].values, data_ch0['base_sig'].values,
                                    guess=[0, 0])
    plt.errorbar(data_ch0['time'].values, data_ch0['base_mean'].values, yerr=data_ch0['base_sig'].values,
                 fmt='.', label='Exposed', markersize=3, capsize=1, linewidth=1, capthick=1)
    plt.plot([data_ch0['time'].values[0], data_ch0['time'].values[-1]],
             linear([data_ch0['time'].values[0], data_ch0['time'].values[-1]], *exposed_fit["popt"]), 'k-', zorder=10)

    data_ch1 = data[data["pmt"] == 1]
    control_fit = fit_straight_line(data_ch1['time'].values, data_ch1['base_mean'].values, data_ch1['base_sig'].values,
                                    guess=[0, 0])
    plt.errorbar(data_ch1['time'].values, data_ch1['base_mean'].values, yerr=data_ch1['base_sig'].values,
                 fmt='.', label='Control', markersize=3, capsize=1, linewidth=1, capthick=1)
    plt.plot([data_ch1['time'].values[0], data_ch1['time'].values[-1]],
             linear([data_ch1['time'].values[0], data_ch1['time'].values[-1]], *control_fit["popt"]), 'k-', zorder=10)

    plt.xlabel('Days from 1% He Onset')
    plt.ylabel('Voltage /mV')
    plt.fill_between([-100, 0], [1000, 1000], alpha=0.1,
                     facecolor='green', label='Atmospheric He')
    plt.fill_between([0, 98], [1000, 1000], alpha=0.1,
                     facecolor='blue', label='1% He')
    plt.fill_between([98, 500], [1000, 1000], alpha=0.1,
                     facecolor='red', label='10% He')
    plt.xlim(-40, 430)
    plt.ylim(978, 981)

    plt.legend(loc='best', bbox_to_anchor=(1.0, 1))
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/PMT_Project/baseline_drift.pdf")
    plt.close(fig1)

    plt.figure(figsize=figsize)
    plt.xlabel('Days from 1% He Onset')
    plt.ylabel('PCT Change /%')
    returns_ch0 = data_ch0['base_mean'].pct_change(1)*100
    returns_ch0.plot(label='Exposed')
    returns_ch1 = data_ch1['base_mean'].pct_change(1)*100
    returns_ch1.plot(label='Control')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/PMT_Project/baseline_pct.pdf")
    plt.close()

    dataframe = pd.concat([pd.DataFrame(data_ch0['base_mean'].values).shift(1),
                           pd.DataFrame(data_ch0['base_mean'].values)], axis=1)
    dataframe.columns = ['t-1', 't+1']
    result0 = dataframe.corr()
    print(result0)
    dataframe = pd.concat([pd.DataFrame(data_ch1['base_mean'].values).shift(1),
                           pd.DataFrame(data_ch1['base_mean'].values)], axis=1)
    dataframe.columns = ['t-1', 't+1']
    result1 = dataframe.corr()
    print(result1)

    plt.figure(figsize=figsize)
    lag_plot(data_ch0['base_mean'], label=r'Exposed $\rho$={:.2f}'.format(result0['t-1'][1]), alpha=0.5)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/PMT_Project/baseline_cor_exposed.pdf")

    plt.figure(figsize=figsize)
    lag_plot(data_ch1['base_mean'], label=r'Control $\rho$={:.2f}'.format(result1['t-1'][1]), alpha=0.5)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig("/Users/williamquinn/Desktop/PMT_Project/baseline_cor_control.pdf")
    '''autocorrelation_plot(data_ch0['base_mean'])
    plt.grid()
    plt.show()'''

    '''# create lagged dataset
    values = pd.DataFrame(data_ch0['base_mean'].values)
    dataframe = pd.concat([values.shift(1), values], axis=1)
    dataframe.columns = ['t-1', 't+1']
    # split into train and test sets
    X = dataframe.values
    train, test = X[1:len(X) - 7], X[len(X) - 7:]
    train_X, train_y = train[:, 0], train[:, 1]
    test_X, test_y = test[:, 0], test[:, 1]

    # persistence model
    def model_persistence(x):
        return x

    # walk-forward validation
    predictions = []
    for x in test_X:
        yhat = model_persistence(x)
        predictions.append(yhat)
    test_score = mean_squared_error(test_y, predictions)
    print('Test MSE: %.3f' % test_score)
    # plot predictions vs expected
    plt.plot(test_y, label='Test')
    plt.plot(predictions, label='Prediction')
    plt.legend(loc='best')
    plt.show()'''

    '''X = data_ch0['base_mean'].values
    train, test = X[1:len(X) - 7], X[len(X) - 7:]
    # train autoregression
    model = AutoReg(train, lags=7)
    model_fit = model.fit()
    print('Coefficients: %s' % model_fit.params)
    # make predictions
    predictions = model_fit.predict(start=len(train), end=len(train) + len(test) - 1, dynamic=False)
    for i in range(len(predictions)):
        print('predicted=%f, expected=%f' % (predictions[i], test[i]))
    rmse = np.sqrt(mean_squared_error(test, predictions))
    print('Test RMSE: %.3f' % rmse)
    # plot results
    plt.plot(test, label='test')
    plt.plot(predictions, label='prediction')
    plt.legend()
    plt.show()'''

    '''data = pd.read_csv("/Users/williamquinn/Desktop/PMT_Project/res_vs_time.csv", header=0,
                       dtype={0: int, 1: int, 2: float, 3: float, 4: float, 5: int, 6: float, 7: float},
                       engine='python')
    data['date'] = data['date'].astype(str)
    data = data[data["pmt"] == 0]
    f = []
    for date in data['date'].values:
        new_date = date[4] + date[5] + '/' + date[2] + date[3] + '/' + '20' + date[0] + date[1]
        f.append(new_date)
    data['date'] = pd.DatetimeIndex(f, dayfirst=True)
    data.set_index('date', inplace=True)
    data = data[~data.index.duplicated()]
    data = data.asfreq('D')
    data["base_mean"] = data["base_mean"].fillna(data["base_mean"].mean())'''

    '''decomposed = seasonal_decompose(data['base_mean'], model='additive')
    decomposed.plot()  # See note below about this
    plt.show()'''

    # data['stationary'] = data['base_mean'].diff()

    '''decomposed = seasonal_decompose(data['stationary'][1:], model='additive')
    decomposed.plot()  # See note below about this
    plt.show()'''

    '''# create train/test datasets
    X = data['stationary'].dropna()
    train_data = X[1:len(X) - 12]
    test_data = X[len(X) - 12:]
    # train the autoregression model
    model = AR(train_data)
    model_fitted = model.fit()
    print('The lag value chose is: %s' % model_fitted.k_ar)
    print('The coefficients of the model are:\n %s' % model_fitted.params)

    # make predictions
    predictions = model_fitted.predict(
        start=len(train_data),
        end=len(train_data) + len(test_data) - 1,
        dynamic=False)

    # create a comparison dataframe
    compare_df = pd.concat(
        [data['stationary'].tail(12),
         predictions], axis=1).rename(
        columns={'stationary': 'actual', 0: 'predicted'})
    # plot the two values
    # compare_df.plot()
    # plt.show()

    r2 = r2_score(data['stationary'].tail(12), predictions)
    print(r2)'''


def main():
    topology = [2, 1]
    pmt_array = PMT_Array(topology, "summary")
    pmt_array.set_pmt_id("GAO607", 0)
    pmt_array.set_pmt_id("GAO612", 1)

    # store_res(pmt_array)
    plot_res()
    plot_base_drift()


if __name__ == '__main__':
    main()
