import numpy as np
from scipy.optimize import curve_fit
from array import array
import ROOT
import pandas as pd


def io_parse_arguments():
    import argparse
    parser = argparse.ArgumentParser(description="Input file names")
    parser.add_argument('-i', required=False, type=str, help='Input data file path')
    parser.add_argument('-o', required=False, type=str, help='Output data file path')
    parser.add_argument('-t', required=False, type=str, help='Output data file path')
    args = parser.parse_args()
    return args


def get_normalisation_factor(vector: list):
    return np.sqrt(np.dot(vector, vector))


def get_n_sfs(popt, perr):
    ns = []
    for p, err in zip(popt, perr):
        p, err = abs(p), abs(err)
        # this will not be in standard form
        if 0.01 < p < 1000:
            if err < 1:
                n = abs(int(np.log10(err))) + 1
            else:
                n = abs(int(np.log10(err)))
        else:
            #if err < 1 and p < 1:
            n = int(np.log10(p)) - int(np.log10(err))
            '''else:
                n = abs(int(np.log10(p))) - abs(int(np.log10(err)))'''
        ns.append(n)
    return ns


def get_run_number(input_data_path: str):
    i = input_data_path.split("/")
    ii = i[-1]
    iii = ii.split("_")
    iv = iii[2]
    v = iv.split(".")
    vi = v[0]
    return str(vi.strip())


def get_voltage(input_data_path: str):
    i = input_data_path.split("/")
    ii = i[-1]
    iii = ii.split("_")
    iv = iii[0]
    v = iv.split("A")
    vi = v[1]
    return str(vi.strip())


def process_date(date_array: np.array):
    output_list = []
    for index, date in enumerate(date_array):
        temp_date = 0
        output_date = 0
        if 190000 < date < 200000:
            temp_date = date - 190000
            if 100 < temp_date < 200:
                temp_date = temp_date - 100
            elif 200 < temp_date < 300:
                output_date += 31
                temp_date = temp_date - 200
            elif 300 < temp_date < 400:
                output_date += 28 + 31
                temp_date = temp_date - 300
            elif 400 < temp_date < 500:
                output_date += 31 + 28 + 31
                temp_date = temp_date - 400
            elif 500 < temp_date < 600:
                output_date += 31 + 28 + 31 + 30
                temp_date = temp_date - 500
            elif 600 < temp_date < 700:
                output_date += 31 + 28 + 31 + 30 + 31
                temp_date = temp_date - 600
            elif 700 < temp_date < 800:
                output_date += 31 + 28 + 31 + 30 + 31 + 30
                temp_date = temp_date - 700
            elif 800 < temp_date < 900:
                output_date += 31 + 28 + 31 + 30 + 31 + 30 + 31
                temp_date = temp_date - 800
            elif 900 < temp_date < 1000:
                output_date += 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31
                temp_date = temp_date - 900
            elif 1000 < temp_date < 1100:
                output_date += 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30
                temp_date = temp_date - 1000
            elif 1100 < temp_date < 1200:
                output_date += 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31
                temp_date = temp_date - 1100
            elif 1200 < temp_date < 1300:
                output_date += 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 + 30
                temp_date = temp_date - 1200
            output_date += temp_date
        elif 200000 < date < 210000:
            temp_date = date - 200000
            output_date += 365
            if 100 < temp_date < 200:
                temp_date = temp_date - 100
            elif 200 < temp_date < 300:
                output_date += 31
                temp_date = temp_date - 200
            elif 300 < temp_date < 400:
                output_date += 29 + 31
                temp_date = temp_date - 300
            elif 400 < temp_date < 500:
                output_date += 31 + 29 + 31
                temp_date = temp_date - 400
            elif 500 < temp_date < 600:
                output_date += 31 + 29 + 31 + 30
                temp_date = temp_date - 500
            elif 600 < temp_date < 700:
                output_date += 31 + 29 + 31 + 30 + 31
                temp_date = temp_date - 600
            elif 700 < temp_date < 800:
                output_date += 31 + 29 + 31 + 30 + 31 + 30
                temp_date = temp_date - 700
            elif 800 < temp_date < 900:
                output_date += 31 + 29 + 31 + 30 + 31 + 30 + 31
                temp_date = temp_date - 800
            elif 900 < temp_date < 1000:
                output_date += 31 + 29 + 31 + 30 + 31 + 30 + 31 + 31
                temp_date = temp_date - 900
            elif 1000 < temp_date < 1100:
                output_date += 31 + 29 + 31 + 30 + 31 + 30 + 31 + 31 + 30
                temp_date = temp_date - 1000
            elif 1100 < temp_date < 1200:
                output_date += 31 + 29 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31
                temp_date = temp_date - 1100
            elif 1200 < temp_date < 1300:
                output_date += 31 + 29 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 + 30
                temp_date = temp_date - 1200
            output_date += temp_date
        elif 210000 < date < 220000:
            temp_date = date - 210000
            output_date += 365 * 2 + 1
            if 100 < temp_date < 200:
                temp_date = temp_date - 100
            elif 200 < temp_date < 300:
                output_date += 31
                temp_date = temp_date - 200
            elif 300 < temp_date < 400:
                output_date += 28 + 31
                temp_date = temp_date - 300
            elif 400 < temp_date < 500:
                output_date += 31 + 28 + 31
                temp_date = temp_date - 400
            elif 500 < temp_date < 600:
                output_date += 31 + 28 + 31 + 30
                temp_date = temp_date - 500
            elif 600 < temp_date < 700:
                output_date += 31 + 28 + 31 + 30 + 31
                temp_date = temp_date - 600
            elif 700 < temp_date < 800:
                output_date += 31 + 28 + 31 + 30 + 31 + 30
                temp_date = temp_date - 700
            elif 800 < temp_date < 900:
                output_date += 31 + 28 + 31 + 30 + 31 + 30 + 31
                temp_date = temp_date - 800
            elif 900 < temp_date < 1000:
                output_date += 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31
                temp_date = temp_date - 900
            elif 1000 < temp_date < 1100:
                output_date += 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30
                temp_date = temp_date - 1000
            elif 1100 < temp_date < 1200:
                output_date += 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31
                temp_date = temp_date - 1100
            elif 1200 < temp_date < 1300:
                output_date += 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 + 30
                temp_date = temp_date - 1200
            output_date += temp_date
        elif 220000 < date < 230000:
            temp_date = date - 220000
            output_date += 365 * 3 + 1
            if 100 < temp_date < 200:
                temp_date = temp_date - 100
            elif 200 < temp_date < 300:
                output_date += 31
                temp_date = temp_date - 200
            elif 300 < temp_date < 400:
                output_date += 28 + 31
                temp_date = temp_date - 300
            elif 400 < temp_date < 500:
                output_date += 31 + 28 + 31
                temp_date = temp_date - 400
            elif 500 < temp_date < 600:
                output_date += 31 + 28 + 31 + 30
                temp_date = temp_date - 500
            elif 600 < temp_date < 700:
                output_date += 31 + 28 + 31 + 30 + 31
                temp_date = temp_date - 600
            elif 700 < temp_date < 800:
                output_date += 31 + 28 + 31 + 30 + 31 + 30
                temp_date = temp_date - 700
            elif 800 < temp_date < 900:
                output_date += 31 + 28 + 31 + 30 + 31 + 30 + 31
                temp_date = temp_date - 800
            elif 900 < temp_date < 1000:
                output_date += 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31
                temp_date = temp_date - 900
            elif 1000 < temp_date < 1100:
                output_date += 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30
                temp_date = temp_date - 1000
            elif 1100 < temp_date < 1200:
                output_date += 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31
                temp_date = temp_date - 1100
            elif 1200 < temp_date < 1300:
                output_date += 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 + 30
                temp_date = temp_date - 1200
            output_date += temp_date
        output_date -= 310
        output_list.append(output_date)

    output_array = np.array(output_list)
    assert output_array.size == len(date_array)
    return output_array


def process_exposure(date_array: np.array):
    output_list = []
    for index, date in enumerate(date_array):
        temp_date = 0
        if date == 191008:
            temp_date = 0 / 1000000
        if date == 191009:
            temp_date = 1 / 1000000
        if date == 191010:
            temp_date = 2 / 1000000
        if date == 191011:
            temp_date = 3 / 1000000
        if date == 191012:
            temp_date = 4 / 1000000
        if date == 191013:
            temp_date = 5 / 1000000
        if date == 191014:
            temp_date = 6 / 1000000
        if date == 191015:
            temp_date = 7 / 1000000
        if date == 191016:
            temp_date = 8 / 1000000
        if date == 191017:
            temp_date = 9 / 1000000
        if date == 191018:
            temp_date = 10 / 1000000
        if date == 191019:
            temp_date = 11 / 1000000
        if date == 191020:
            temp_date = 12 / 1000000
        if date == 191021:
            temp_date = 13 / 1000000
        if date == 191022:
            temp_date = 14 / 1000000
        if date == 191023:
            temp_date = 15 / 1000000
        if date == 191024:
            temp_date = 16 / 1000000
        if date == 191025:
            temp_date = 17 / 1000000
        if date == 191026:
            temp_date = 18 / 1000000
        if date == 191027:
            temp_date = 19 / 1000000
        if date == 191028:
            temp_date = 20 / 1000000
        if date == 191029:
            temp_date = 21 / 1000000
        if date == 191030:
            temp_date = 22 / 1000000
        if date == 191031:
            temp_date = 23 / 1000000
        if date == 191101:
            temp_date = 24 / 1000000
        if date == 191102:
            temp_date = 25 / 1000000
        if date == 191103:
            temp_date = 26 / 1000000
        if date == 191104:
            temp_date = 27 / 1000000
        if date == 191105:
            temp_date = 28 / 1000000
        if date == 191106:
            temp_date = 29 / 1000000
        if date == 191107:
            temp_date = 1 / 100
        if date == 191108:
            temp_date = 2 / 100
        if date == 191109:
            temp_date = 3 / 100
        if date == 191110:
            temp_date = 4 / 100
        if date == 191111:
            temp_date = 5 / 100
        if date == 191112:
            temp_date = 6 / 100
        if date == 191113:
            temp_date = 7 / 100
        if date == 191114:
            temp_date = 8 / 100
        if date == 191115:
            temp_date = 9 / 100
        if date == 191116:
            temp_date = 10 / 100
        if date == 191117:
            temp_date = 11 / 100
        if date == 191118:
            temp_date = 12 / 100
        if date == 191119:
            temp_date = 13 / 100
        if date == 191120:
            temp_date = 14 / 100
        if date == 191121:
            temp_date = 15 / 100
        if date == 191122:
            temp_date = 16 / 100
        if date == 191123:
            temp_date = 17 / 100
        if date == 191124:
            temp_date = 18 / 100
        if date == 191125:
            temp_date = 19 / 100
        if date == 191126:
            temp_date = 20 / 100
        if date == 191127:
            temp_date = 21 / 100
        if date == 191128:
            temp_date = 22 / 100
        if date == 191129:
            temp_date = 23 / 100
        if date == 191130:
            temp_date = 24 / 100
        if date == 191201:
            temp_date = 25 / 100
        if date == 191202:
            temp_date = 26 / 100
        if date == 191203:
            temp_date = 27 / 100
        if date == 191204:
            temp_date = 28 / 100
        if date == 191205:
            temp_date = 29 / 100
        if date == 191206:
            temp_date = 30 / 100
        if date == 191207:
            temp_date = 31 / 100
        if date == 191208:
            temp_date = 32 / 100
        if date == 191209:
            temp_date = 33 / 100
        if date == 191210:
            temp_date = 34 / 100
        if date == 191211:
            temp_date = 35 / 100
        if date == 191212:
            temp_date = 36 / 100
        if date == 191213:
            temp_date = 37 / 100
        if date == 191214:
            temp_date = 38 / 100
        if date == 191215:
            temp_date = 39 / 100
        if date == 191216:
            temp_date = 40 / 100
        if date == 191217:
            temp_date = 41 / 100
        if date == 191218:
            temp_date = 42 / 100
        if date == 191219:
            temp_date = 43 / 100
        if date == 191220:
            temp_date = 44 / 100
        if date == 191221:
            temp_date = 45 / 100
        if date == 191222:
            temp_date = 46 / 100
        if date == 191223:
            temp_date = 47 / 100
        if date == 191224:
            temp_date = 48 / 100
        if date == 191225:
            temp_date = 49 / 100
        if date == 191226:
            temp_date = 50 / 100
        if date == 191227:
            temp_date = 51 / 100
        if date == 191228:
            temp_date = 52 / 100
        if date == 191229:
            temp_date = 53 / 100
        if date == 191230:
            temp_date = 54 / 100
        if date == 191231:
            temp_date = 55 / 100
        if date == 200101:
            temp_date = 56 / 100
        if date == 200102:
            temp_date = 57 / 100
        if date == 200103:
            temp_date = 58 / 100
        if date == 200104:
            temp_date = 59 / 100
        if date == 200105:
            temp_date = 60 / 100
        if date == 200106:
            temp_date = 61 / 100
        if date == 200107:
            temp_date = 62 / 100
        if date == 200108:
            temp_date = 63 / 100
        if date == 200109:
            temp_date = 64 / 100
        if date == 200110:
            temp_date = 65 / 100
        if date == 200111:
            temp_date = 66 / 100
        if date == 200112:
            temp_date = 67 / 100
        if date == 200113:
            temp_date = 68 / 100
        if date == 200114:
            temp_date = 69 / 100
        if date == 200115:
            temp_date = 70 / 100
        if date == 200116:
            temp_date = 71 / 100
        if date == 200117:
            temp_date = 72 / 100
        if date == 200118:
            temp_date = 73 / 100
        if date == 200119:
            temp_date = 74 / 100
        if date == 200120:
            temp_date = 75 / 100
        if date == 200121:
            temp_date = 76 / 100
        if date == 200122:
            temp_date = 77 / 100
        if date == 200123:
            temp_date = 78 / 100
        if date == 200124:
            temp_date = 79 / 100
        if date == 200125:
            temp_date = 80 / 100
        if date == 200126:
            temp_date = 81 / 100
        if date == 200127:
            temp_date = 82 / 100
        if date == 200128:
            temp_date = 83 / 100
        if date == 200129:
            temp_date = 84 / 100
        if date == 200130:
            temp_date = 85 / 100
        if date == 200131:
            temp_date = 86 / 100
        if date == 200201:
            temp_date = 87 / 100
        if date == 200202:
            temp_date = 88 / 100
        if date == 200203:
            temp_date = 89 / 100
        if date == 200204:
            temp_date = 90 / 100
        if date == 200205:
            temp_date = 91 / 100
        if date == 200206:
            temp_date = 92 / 100
        if date == 200207:
            temp_date = 93 / 100
        if date == 200208:
            temp_date = 94 / 100
        if date == 200209:
            temp_date = 95 / 100
        if date == 200210:
            temp_date = 96 / 100
        if date == 200211:
            temp_date = 97 / 100
        if date == 200212:
            temp_date = 1 / 10 + 97 / 100
        if date == 200213:
            temp_date = 2 / 10 + 97 / 100
        if date == 200214:
            temp_date = 3 / 10 + 97 / 100
        if date == 200215:
            temp_date = 4 / 10 + 97 / 100
        if date == 200216:
            temp_date = 5 / 10 + 97 / 100
        if date == 200217:
            temp_date = 6 / 10 + 97 / 100
        if date == 200218:
            temp_date = 7 / 10 + 97 / 100
        if date == 200219:
            temp_date = 8 / 10 + 97 / 100
        if date == 200220:
            temp_date = 9 / 10 + 97 / 100
        if date == 200220:
            temp_date = 10 / 10 + 97 / 100
        if date == 200221:
            temp_date = 11 / 10 + 97 / 100
        if date == 200222:
            temp_date = 12 / 10 + 97 / 100
        if date == 200223:
            temp_date = 13 / 10 + 97 / 100
        if date == 200224:
            temp_date = 14 / 10 + 97 / 100
        if date == 200225:
            temp_date = 15 / 10 + 97 / 100
        if date == 200226:
            temp_date = 16 / 10 + 97 / 100
        if date == 200227:
            temp_date = 17 / 10 + 97 / 100
        if date == 200228:
            temp_date = 18 / 10 + 97 / 100
        if date == 200229:
            temp_date = 19 / 10 + 97 / 100
        if date == 200301:
            temp_date = 20 / 10 + 97 / 100
        if date == 200302:
            temp_date = 21 / 10 + 97 / 100
        if date == 200303:
            temp_date = 22 / 10 + 97 / 100
        if date == 200304:
            temp_date = 23 / 10 + 97 / 100
        if date == 200305:
            temp_date = 24 / 10 + 97 / 100
        if date == 200306:
            temp_date = 25 / 10 + 97 / 100
        if date == 200307:
            temp_date = 26 / 10 + 97 / 100
        if date == 200308:
            temp_date = 27 / 10 + 97 / 100
        if date == 200309:
            temp_date = 28 / 10 + 97 / 100
        if date == 200310:
            temp_date = 29 / 10 + 97 / 100
        if date == 200311:
            temp_date = 30 / 10 + 97 / 100
        if date == 200312:
            temp_date = 31 / 10 + 97 / 100
        if date == 200313:
            temp_date = 32 / 10 + 97 / 100
        if date == 200314:
            temp_date = 33 / 10 + 97 / 100
        if date == 200315:
            temp_date = 34 / 10 + 97 / 100
        if date == 200316:
            temp_date = 35 / 10 + 97 / 100
        if date == 200317:
            temp_date = 36 / 10 + 97 / 100
        if date == 200318:
            temp_date = 37 / 10 + 97 / 100
        if date == 200319:
            temp_date = 38 / 10 + 97 / 100
        if date == 200320:
            temp_date = 39 / 10 + 97 / 100
        if date == 200321:
            temp_date = 40 / 10 + 97 / 100
        if date == 200322:
            temp_date = 41 / 10 + 97 / 100
        if date == 200323:
            temp_date = 42 / 10 + 97 / 100
        if date == 200324:
            temp_date = 43 / 10 + 97 / 100
        if date == 200325:
            temp_date = 44 / 10 + 97 / 100
        if date == 200326:
            temp_date = 45 / 10 + 97 / 100
        if date == 200327:
            temp_date = 46 / 10 + 97 / 100
        if date == 200328:
            temp_date = 47 / 10 + 97 / 100
        if date == 200329:
            temp_date = 48 / 10 + 97 / 100
        if date == 200330:
            temp_date = 49 / 10 + 97 / 100
        if date == 200331:
            temp_date = 50 / 10 + 97 / 100
        if date == 200401:
            temp_date = 51 / 10 + 97 / 100
        if date == 200402:
            temp_date = 52 / 10 + 97 / 100
        if date == 200403:
            temp_date = 53 / 10 + 97 / 100
        if date == 200404:
            temp_date = 54 / 10 + 97 / 100
        if date == 200405:
            temp_date = 55 / 10 + 97 / 100
        if date == 200406:
            temp_date = 56 / 10 + 97 / 100
        if date == 200407:
            temp_date = 57 / 10 + 97 / 100
        if date == 200408:
            temp_date = 58 / 10 + 97 / 100
        if date == 200409:
            temp_date = 59 / 10 + 97 / 100
        if date == 200410:
            temp_date = 60 / 10 + 97 / 100
        if date == 200411:
            temp_date = 61 / 10 + 97 / 100
        if date == 200412:
            temp_date = 62 / 10 + 97 / 100
        if date == 200413:
            temp_date = 63 / 10 + 97 / 100
        if date == 200414:
            temp_date = 64 / 10 + 97 / 100
        if date == 200415:
            temp_date = 65 / 10 + 97 / 100
        if date == 200416:
            temp_date = 66 / 10 + 97 / 100
        if date == 200417:
            temp_date = 67 / 10 + 97 / 100
        if date == 200418:
            temp_date = 68 / 10 + 97 / 100
        if date == 200419:
            temp_date = 69 / 10 + 97 / 100
        if date == 200420:
            temp_date = 70 / 10 + 97 / 100
        if date == 200421:
            temp_date = 71 / 10 + 97 / 100
        if date == 200422:
            temp_date = 72 / 10 + 97 / 100
        if date == 200423:
            temp_date = 73 / 10 + 97 / 100
        if date == 200424:
            temp_date = 74 / 10 + 97 / 100
        if date == 200425:
            temp_date = 75 / 10 + 97 / 100
        if date == 200426:
            temp_date = 76 / 10 + 97 / 100
        if date == 200427:
            temp_date = 77 / 10 + 97 / 100
        if date == 200428:
            temp_date = 78 / 10 + 97 / 100
        if date == 200429:
            temp_date = 79 / 10 + 97 / 100
        if date == 200430:
            temp_date = 80 / 10 + 97 / 100
        if date == 200501:
            temp_date = 81 / 10 + 97 / 100
        if date == 200502:
            temp_date = 82 / 10 + 97 / 100
        if date == 200503:
            temp_date = 83 / 10 + 97 / 100
        if date == 200504:
            temp_date = 84 / 10 + 97 / 100
        if date == 200505:
            temp_date = 85 / 10 + 97 / 100
        if date == 200506:
            temp_date = 86 / 10 + 97 / 100
        if date == 200507:
            temp_date = 87 / 10 + 97 / 100
        if date == 200508:
            temp_date = 88 / 10 + 97 / 100
        if date == 200509:
            temp_date = 89 / 10 + 97 / 100
        if date == 200510:
            temp_date = 90 / 10 + 97 / 100
        if date == 200511:
            temp_date = 91 / 10 + 97 / 100
        if date == 200512:
            temp_date = 92 / 10 + 97 / 100
        if date == 200513:
            temp_date = 93 / 10 + 97 / 100
        if date == 200514:
            temp_date = 94 / 10 + 97 / 100
        if date == 200515:
            temp_date = 95 / 10 + 97 / 100
        if date == 200516:
            temp_date = 96 / 10 + 97 / 100
        if date == 200517:
            temp_date = 97 / 10 + 97 / 100
        if date == 200518:
            temp_date = 98 / 10 + 97 / 100
        if date == 200519:
            temp_date = 99 / 10 + 97 / 100
        if date == 200520:
            temp_date = 100 / 10 + 97 / 100
        if date == 200521:
            temp_date = 101 / 10 + 97 / 100
        if date == 200522:
            temp_date = 102 / 10 + 97 / 100
        if date == 200523:
            temp_date = 103 / 10 + 97 / 100
        if date == 200524:
            temp_date = 104 / 10 + 97 / 100
        if date == 200525:
            temp_date = 105 / 10 + 97 / 100
        if date == 200526:
            temp_date = 106 / 10 + 97 / 100
        if date == 200527:
            temp_date = 107 / 10 + 97 / 100
        if date == 200528:
            temp_date = 108 / 10 + 97 / 100
        if date == 200529:
            temp_date = 109 / 10 + 97 / 100
        if date == 200530:
            temp_date = 110 / 10 + 97 / 100
        if date == 200531:
            temp_date = 111 / 10 + 97 / 100
        if date == 200601:
            temp_date = 112 / 10 + 97 / 100
        if date == 200602:
            temp_date = 113 / 10 + 97 / 100
        if date == 200603:
            temp_date = 114 / 10 + 97 / 100
        if date == 200604:
            temp_date = 115 / 10 + 97 / 100
        if date == 200605:
            temp_date = 116 / 10 + 97 / 100
        if date == 200606:
            temp_date = 117 / 10 + 97 / 100
        if date == 200607:
            temp_date = 118 / 10 + 97 / 100
        if date == 200608:
            temp_date = 119 / 10 + 97 / 100
        if date == 200609:
            temp_date = 120 / 10 + 97 / 100
        if date == 200610:
            temp_date = 121 / 10 + 97 / 100
        if date == 200611:
            temp_date = 122 / 10 + 97 / 100
        if date == 200612:
            temp_date = 123 / 10 + 97 / 100
        if date == 200613:
            temp_date = 124 / 10 + 97 / 100
        if date == 200614:
            temp_date = 125 / 10 + 97 / 100
        if date == 200615:
            temp_date = 126 / 10 + 97 / 100
        if date == 200616:
            temp_date = 127 / 10 + 97 / 100
        if date == 200617:
            temp_date = 128 / 10 + 97 / 100
        if date == 200618:
            temp_date = 129 / 10 + 97 / 100
        if date == 200619:
            temp_date = 130 / 10 + 97 / 100
        if date == 200620:
            temp_date = 131 / 10 + 97 / 100
        if date == 200621:
            temp_date = 132 / 10 + 97 / 100
        if date == 200622:
            temp_date = 133 / 10 + 97 / 100
        if date == 200623:
            temp_date = 134 / 10 + 97 / 100
        if date == 200624:
            temp_date = 135 / 10 + 97 / 100
        if date == 200625:
            temp_date = 136 / 10 + 97 / 100
        if date == 200626:
            temp_date = 137 / 10 + 97 / 100
        if date == 200627:
            temp_date = 138 / 10 + 97 / 100
        if date == 200628:
            temp_date = 139 / 10 + 97 / 100
        if date == 200629:
            temp_date = 140 / 10 + 97 / 100
        if date == 200630:
            temp_date = 141 / 10 + 97 / 100
        if date == 200701:
            temp_date = 142 / 10 + 97 / 100
        if date == 200702:
            temp_date = 143 / 10 + 97 / 100
        if date == 200703:
            temp_date = 144 / 10 + 97 / 100
        if date == 200704:
            temp_date = 145 / 10 + 97 / 100
        if date == 200705:
            temp_date = 146 / 10 + 97 / 100
        if date == 200706:
            temp_date = 147 / 10 + 97 / 100
        if date == 200707:
            temp_date = 148 / 10 + 97 / 100
        if date == 200708:
            temp_date = 149 / 10 + 97 / 100
        if date == 200709:
            temp_date = 150 / 10 + 97 / 100
        if date == 200710:
            temp_date = 151 / 10 + 97 / 100
        if date == 200711:
            temp_date = 152 / 10 + 97 / 100
        if date == 200712:
            temp_date = 153 / 10 + 97 / 100
        if date == 200713:
            temp_date = 154 / 10 + 97 / 100
        if date == 200714:
            temp_date = 155 / 10 + 97 / 100
        if date == 200715:
            temp_date = 156 / 10 + 97 / 100
        if date == 200716:
            temp_date = 157 / 10 + 97 / 100
        if date == 200717:
            temp_date = 158 / 10 + 97 / 100
        if date == 200718:
            temp_date = 159 / 10 + 97 / 100
        if date == 200719:
            temp_date = 160 / 10 + 97 / 100
        if date == 200720:
            temp_date = 161 / 10 + 97 / 100
        if date == 200721:
            temp_date = 162 / 10 + 97 / 100
        if date == 200722:
            temp_date = 163 / 10 + 97 / 100
        if date == 200723:
            temp_date = 164 / 10 + 97 / 100
        if date == 200724:
            temp_date = 165 / 10 + 97 / 100
        if date == 200725:
            temp_date = 166 / 10 + 97 / 100
        if date == 200726:
            temp_date = 167 / 10 + 97 / 100
        if date == 200727:
            temp_date = 168 / 10 + 97 / 100
        if date == 200728:
            temp_date = 169 / 10 + 97 / 100
        if date == 200729:
            temp_date = 170 / 10 + 97 / 100
        if date == 200730:
            temp_date = 171 / 10 + 97 / 100
        if date == 200731:
            temp_date = 171 / 10 + 97 / 100
        if date == 200801:
            temp_date = 172 / 10 + 97 / 100
        if date == 200802:
            temp_date = 173 / 10 + 97 / 100
        if date == 200803:
            temp_date = 174 / 10 + 97 / 100
        if date == 200804:
            temp_date = 175 / 10 + 97 / 100
        if date == 200805:
            temp_date = 176 / 10 + 97 / 100
        if date == 200806:
            temp_date = 177 / 10 + 97 / 100
        if date == 200807:
            temp_date = 178 / 10 + 97 / 100
        if date == 200808:
            temp_date = 179 / 10 + 97 / 100
        if date == 200809:
            temp_date = 180 / 10 + 97 / 100
        if date == 200810:
            temp_date = 181 / 10 + 97 / 100
        if date == 200811:
            temp_date = 182 / 10 + 97 / 100
        if date == 200812:
            temp_date = 183 / 10 + 97 / 100
        if date == 200813:
            temp_date = 184 / 10 + 97 / 100
        if date == 200814:
            temp_date = 185 / 10 + 97 / 100
        if date == 200815:
            temp_date = 186 / 10 + 97 / 100
        if date == 200816:
            temp_date = 187 / 10 + 97 / 100
        if date == 200817:
            temp_date = 188 / 10 + 97 / 100
        if date == 200818:
            temp_date = 189 / 10 + 97 / 100
        if date == 200819:
            temp_date = 190 / 10 + 97 / 100
        if date == 200820:
            temp_date = 191 / 10 + 97 / 100
        if date == 200821:
            temp_date = 192 / 10 + 97 / 100
        if date == 200822:
            temp_date = 193 / 10 + 97 / 100
        if date == 200823:
            temp_date = 194 / 10 + 97 / 100
        if date == 200824:
            temp_date = 195 / 10 + 97 / 100
        if date == 200825:
            temp_date = 196 / 10 + 97 / 100
        if date == 200826:
            temp_date = 197 / 10 + 97 / 100
        if date == 200827:
            temp_date = 198 / 10 + 97 / 100
        if date == 200828:
            temp_date = 199 / 10 + 97 / 100
        if date == 200829:
            temp_date = 200 / 10 + 97 / 100
        if date == 200830:
            temp_date = 201 / 10 + 97 / 100
        if date == 200831:
            temp_date = 202 / 10 + 97 / 100
        if date == 200901:
            temp_date = 203 / 10 + 97 / 100
        if date == 200902:
            temp_date = 204 / 10 + 97 / 100
        if date == 200903:
            temp_date = 205 / 10 + 97 / 100
        if date == 200904:
            temp_date = 206 / 10 + 97 / 100
        if date == 200905:
            temp_date = 207 / 10 + 97 / 100
        if date == 200906:
            temp_date = 208 / 10 + 97 / 100
        if date == 200907:
            temp_date = 209 / 10 + 97 / 100
        if date == 200908:
            temp_date = 210 / 10 + 97 / 100
        if date == 200909:
            temp_date = 211 / 10 + 97 / 100
        if date == 200910:
            temp_date = 212 / 10 + 97 / 100
        if date == 200911:
            temp_date = 213 / 10 + 97 / 100
        if date == 200912:
            temp_date = 214 / 10 + 97 / 100
        if date == 200913:
            temp_date = 215 / 10 + 97 / 100
        if date == 200914:
            temp_date = 216 / 10 + 97 / 100
        if date == 200915:
            temp_date = 217 / 10 + 97 / 100
        if date == 200916:
            temp_date = 218 / 10 + 97 / 100
        if date == 200917:
            temp_date = 219 / 10 + 97 / 100
        if date == 200918:
            temp_date = 220 / 10 + 97 / 100
        if date == 200919:
            temp_date = 221 / 10 + 97 / 100
        if date == 200920:
            temp_date = 222 / 10 + 97 / 100
        if date == 200921:
            temp_date = 223 / 10 + 97 / 100
        if date == 200922:
            temp_date = 224 / 10 + 97 / 100
        if date == 200923:
            temp_date = 225 / 10 + 97 / 100
        if date == 200924:
            temp_date = 226 / 10 + 97 / 100
        if date == 200925:
            temp_date = 227 / 10 + 97 / 100
        if date == 200926:
            temp_date = 228 / 10 + 97 / 100
        if date == 200927:
            temp_date = 229 / 10 + 97 / 100
        if date == 200928:
            temp_date = 230 / 10 + 97 / 100
        if date == 200929:
            temp_date = 231 / 10 + 97 / 100
        if date == 200930:
            temp_date = 232 / 10 + 97 / 100
        if date == 201001:
            temp_date = 233 / 10 + 97 / 100
        if date == 201002:
            temp_date = 234 / 10 + 97 / 100
        if date == 201003:
            temp_date = 235 / 10 + 97 / 100
        if date == 201004:
            temp_date = 236 / 10 + 97 / 100
        if date == 201005:
            temp_date = 237 / 10 + 97 / 100
        if date == 201006:
            temp_date = 238 / 10 + 97 / 100
        if date == 201007:
            temp_date = 239 / 10 + 97 / 100
        if date == 201008:
            temp_date = 240 / 10 + 97 / 100
        if date == 201009:
            temp_date = 241 / 10 + 97 / 100
        if date == 201010:
            temp_date = 242 / 10 + 97 / 100
        if date == 201011:
            temp_date = 243 / 10 + 97 / 100
        if date == 201012:
            temp_date = 244 / 10 + 97 / 100
        if date == 201013:
            temp_date = 245 / 10 + 97 / 100
        if date == 201014:
            temp_date = 246 / 10 + 97 / 100
        if date == 201015:
            temp_date = 247 / 10 + 97 / 100
        if date == 201016:
            temp_date = 248 / 10 + 97 / 100
        if date == 201017:
            temp_date = 249 / 10 + 97 / 100
        if date == 201018:
            temp_date = 250 / 10 + 97 / 100
        if date == 201019:
            temp_date = 251 / 10 + 97 / 100
        if date == 201020:
            temp_date = 252 / 10 + 97 / 100
        if date == 201021:
            temp_date = 253 / 10 + 97 / 100
        if date == 201022:
            temp_date = 254 / 10 + 97 / 100
        if date == 201023:
            temp_date = 255 / 10 + 97 / 100
        if date == 201024:
            temp_date = 256 / 10 + 97 / 100
        if date == 201025:
            temp_date = 257 / 10 + 97 / 100
        if date == 201026:
            temp_date = 258 / 10 + 97 / 100
        if date == 201027:
            temp_date = 259 / 10 + 97 / 100
        if date == 201028:
            temp_date = 260 / 10 + 97 / 100
        if date == 201029:
            temp_date = 261 / 10 + 97 / 100
        if date == 201030:
            temp_date = 262 / 10 + 97 / 100
        if date == 201031:
            temp_date = 263 / 10 + 97 / 100
        if date == 201101:
            temp_date = 264 / 10 + 97 / 100
        if date == 201102:
            temp_date = 265 / 10 + 97 / 100
        if date == 201103:
            temp_date = 266 / 10 + 97 / 100
        if date == 201104:
            temp_date = 267 / 10 + 97 / 100
        if date == 201105:
            temp_date = 268 / 10 + 97 / 100
        if date == 201106:
            temp_date = 269 / 10 + 97 / 100
        if date == 201107:
            temp_date = 270 / 10 + 97 / 100
        if date == 201108:
            temp_date = 271 / 10 + 97 / 100
        if date == 201109:
            temp_date = 272 / 10 + 97 / 100
        if date == 201110:
            temp_date = 273 / 10 + 97 / 100
        if date == 201111:
            temp_date = 274 / 10 + 97 / 100
        if date == 201112:
            temp_date = 275 / 10 + 97 / 100
        if date == 201113:
            temp_date = 276 / 10 + 97 / 100
        if date == 201114:
            temp_date = 277 / 10 + 97 / 100
        if date == 201115:
            temp_date = 278 / 10 + 97 / 100
        if date == 201116:
            temp_date = 279 / 10 + 97 / 100
        if date == 201117:
            temp_date = 280 / 10 + 97 / 100
        if date == 201118:
            temp_date = 281 / 10 + 97 / 100
        if date == 201119:
            temp_date = 282 / 10 + 97 / 100
        if date == 201120:
            temp_date = 283 / 10 + 97 / 100
        if date == 201121:
            temp_date = 284 / 10 + 97 / 100
        if date == 201122:
            temp_date = 285 / 10 + 97 / 100
        if date == 201123:
            temp_date = 286 / 10 + 97 / 100
        if date == 201124:
            temp_date = 287 / 10 + 97 / 100
        if date == 201125:
            temp_date = 288 / 10 + 97 / 100
        if date == 201126:
            temp_date = 289 / 10 + 97 / 100
        if date == 201127:
            temp_date = 290 / 10 + 97 / 100
        if date == 201128:
            temp_date = 291 / 10 + 97 / 100
        if date == 201129:
            temp_date = 292 / 10 + 97 / 100
        if date == 201130:
            temp_date = 293 / 10 + 97 / 100

        output_list.append(temp_date)

    output_array = np.array(output_list)
    assert output_array.size == date_array.size
    return output_array


def chi2(y_obs, y_err, y_exp, n_par):
    chi2 = 0
    ndof = len(y_obs) - n_par - 1
    for i in range(len(y_exp)):
        chi2 += ((y_exp[i] - y_obs[i]) / y_err[i]) ** 2
    chi2 = chi2
    return chi2, ndof


def gaussian(x, mean, sigma, amplitude, height):
    y = []
    for i in range(len(x)):
        y.append(amplitude * np.exp((-0.5 * ((x[i] - mean) / sigma) ** 2)) + height)
    return y


def gaussian_noh(x, mean, sigma, amplitude):
    y = []
    for i in range(len(x)):
        y.append(amplitude * np.exp((-0.5 * ((x[i] - mean) / sigma) ** 2)))
    return y


def linear(x, m, c):
    y = []
    for i in range(len(x)):
        y.append(m * x[i] + c)
    return y


def om_type(omnum: int):
    string = ''
    index = 0
    if 0 <= omnum < 260:
        row = omnum % 13
        if row == 12 or row == 0:
            string = '5inchMW'
            index = 0
        else:
            string = '8inchMW'
            index = 1
    elif omnum < 520:
        omnum = omnum - 260
        row = omnum % 13
        if row == 12 or row == 0:
            string = '5inchMW'
            index = 0
        else:
            string = '8inchMW'
            index = 1
    elif omnum < 584:
        string = '5inchXW'
        index = 2
    elif omnum < 648:
        string = '5inchXW'
        index = 2
    elif omnum < 680:
        string = '5inchGV'
        index = 3
    elif omnum < 712:
        string = '5inchGV'
        index = 3
    return index, string


def get_om_num(typ, side, wall, col, row):
    if typ == 'M':
        om = row + col*13 + side*260
    elif typ == 'X':
        om = 520 + side*64 + wall*32 + col*16 + row
    elif typ == 'G':
        om = 520 + 128 + side*32 + wall*16 + col
    else:
        raise ValueError("Type of OM needs to be M, X or G")
    return om


def om_id_string(omnum: int):
    string = ''
    if 0 <= omnum < 260:
        col = omnum // 13
        row = omnum % 13
        string = f'M:0.{col}.{row}'
    elif omnum < 520:
        omnum = omnum - 260
        col = omnum // 13
        row = omnum % 13
        string = f'M:1.{col}.{row}'
    elif omnum < 584:
        omnum = omnum - 520
        wall = omnum // 32
        col = omnum % 32 // 16
        row = omnum % 32 % 16
        string = f'X:0.{wall}.{col}.{row}'
    elif omnum < 648:
        omnum = omnum - 520 - 64
        wall = omnum // 32
        col = omnum % 32 // 16
        row = omnum % 32 % 16
        string = f'X:1.{wall}.{col}.{row}'
    elif omnum < 680:
        omnum = omnum - 520 - 128
        wall = omnum % 16
        col = omnum // 16
        string = f'G:0.{col}.{wall}'
    elif omnum < 712:
        omnum = omnum - 520 - 128 - 32
        wall = omnum % 16
        col = omnum // 16
        string = f'G:1.{col}.{wall}'
    return string


def om_id(omnum: int):
    id_ = {"side": None, "wall": None, "col": None, 'row': None}
    if 0 <= omnum < 260:
        id_['side'] = 0
        id_['col'] = omnum // 13
        id_['row'] = omnum % 13
    elif omnum < 520:
        omnum = omnum - 260
        id_['side'] = 1
        id_['col'] = omnum // 13
        id_['row'] = omnum % 13
    elif omnum < 584:
        omnum = omnum - 520
        id_['side'] = 0
        id_['wall'] = omnum // 32
        id_['col'] = omnum % 32 // 16
        id_['row'] = omnum % 32 % 16
    elif omnum < 648:
        omnum = omnum - 520 - 64
        id_['side'] = 1
        id_['wall'] = omnum // 32
        id_['col'] = omnum % 32 // 16
        id_['row'] = omnum % 32 % 16
    elif omnum < 680:
        omnum = omnum - 520 - 128
        id_['side'] = 0
        id_['wall'] = omnum % 16
        id_['col'] = omnum // 16
    elif omnum < 712:
        omnum = omnum - 520 - 128 - 32
        id_['side'] = 1
        id_['wall'] = omnum % 16
        id_['col'] = omnum // 16
    return id_


def get_pmt_type(omnum: id):
    pmt_type = 0
    if 0 <= omnum < 260:
        col = omnum // 13
        row = omnum % 13
        if row == 0 or row == 12:
            pmt_type = 5
        else:
            pmt_type = 8
    elif omnum < 520:
        omnum = omnum - 260
        col = omnum // 13
        row = omnum % 13
        if row == 0 or row == 12:
            pmt_type = 5
        else:
            pmt_type = 8
    else:
        pmt_type = 5

    return pmt_type


def get_cell_id(slot, channel):
    if slot == 0:
        if 0 <= channel < 9:
            side = 1
            row = slot * 2 + 1
            layer = channel
        elif 9 <= channel < 18:
            side = 1
            row = slot * 2
            layer = channel - 9
        elif 18 <= channel < 27:
            side = 0
            row = slot * 2
            layer = channel - 18
        else:
            side = 0
            row = slot * 2 + 1
            layer = channel - 27
    else:
        if 0 <= channel < 9:
            side = 0
            row = slot * 2
            layer = channel
        elif 9 <= channel < 18:
            side = 0
            row = slot * 2 + 1
            layer = channel - 9
        elif 18 <= channel < 27:
            side = 1
            row = slot * 2 + 1
            layer = channel - 18
        else:
            side = 1
            row = slot * 2
            layer = channel - 27
    return '{}.{}.{}'.format(side, row, layer)


def cell_id(cellnum):
    cell_side = cellnum // (9 * 113)
    cell_row = cellnum % (9 * 113) // 9
    cell_layer = cellnum % (9 * 113) % 9

    return cell_side, cell_row, cell_layer


def do_bi_1400(hist, mean):
    lower = mean - 10
    higher = mean + 45

    fit = ROOT.TF1("fit",
                   "[0]*"
                   "(7.08*TMath::Gaus(x,[1],[2]) "
                   " + 1.84*TMath::Gaus(x,[1]*(1 + 72.144/975.651),[2]*1.036) "
                   " + 0.44*TMath::Gaus(x,[1]*(1 + 84.154/975.651),[2]*1.042)) ",
                   lower, higher)

    fit.SetParNames("A", "mu", "sigma")

    fit.SetParLimits(0, 0, 400)
    fit.SetParLimits(1, lower, higher)
    fit.SetParLimits(2, 0.8, 100)
    fit.SetParameters(319, (higher + lower) / 2, 1.09)

    hist.Fit("fit", "0Q", "", lower, higher)

    mu = fit.GetParameter(1)
    mu_err = fit.GetParError(1)
    sig = fit.GetParameter(2)
    sig_err = fit.GetParError(2)
    A = fit.GetParameter(0)
    A_err = fit.GetParError(0)
    chi = fit.GetChisquare() / fit.GetNDF()
    ndof = fit.GetNDF()
    del fit
    return [[mu, mu_err], [sig, sig_err], [A, A_err], [chi*ndof, ndof]]


def do_bi_1000(hist, mean):
    lower = mean - 1
    higher = mean + 6

    fit = ROOT.TF1("fit",
                   "[0]*"
                   "(7.08*TMath::Gaus(x,[1],[2]) "
                   " + 1.84*TMath::Gaus(x,[1]*(1 + 72.144/975.651),[2]*1.036) "
                   " + 0.44*TMath::Gaus(x,[1]*(1 + 84.154/975.651),[2]*1.042)) ",
                   lower, higher)

    fit.SetParNames("A", "mu", "sigma")

    fit.SetParLimits(0, 0, 400)
    fit.SetParLimits(1, lower, higher)
    fit.SetParLimits(2, 0.8, 10)
    fit.SetParameters(319, (higher + lower) / 2, 1.09)

    hist.Fit("fit", "0Q", "", lower, higher)

    mu = fit.GetParameter(1)
    mu_err = fit.GetParError(1)
    sig = fit.GetParameter(2)
    sig_err = fit.GetParError(2)
    A = fit.GetParameter(0)
    A_err = fit.GetParError(0)
    ndof = fit.GetNDF()
    if ndof == 0:
        chi = -1
    else:
        chi = fit.GetChisquare() / fit.GetNDF()
    del fit
    return [[mu, mu_err], [sig, sig_err], [A, A_err], [chi*ndof, ndof]]


def get_sat_bi_mean(hist):
    max_bin = 0
    max_val = 0
    for i_bin in range(50, hist.GetNbinsX() + 1):
        if hist.GetBinContent(i_bin) > max_val:
            max_val = hist.GetBinContent(i_bin)
            max_bin = hist.GetBinLowEdge(i_bin)
    return max_bin


def my_matrix():
    df = pd.read_csv("/Users/williamquinn/Desktop/PMT_Project/HeMatrix.txt",
                     index_col=0)
    df = df / 500
    df.loc["i0", "0"] = 1
    M = df.values

    return M


def quadratic(x, q, m, c):
    return q * x * x + m * x + c


def adjust_pars():
    averages = []
    stds = []
    M = my_matrix()
    nums = np.array([i for i in range(len(M[0]))])
    for i_row in range(len(M)):
        weights = []
        for i_col in range(len(M[i_row])):
            weights.append(M[i_row][i_col])
        weights = np.array(weights)
        average = np.average(nums, weights=weights)

        std = np.sqrt(np.sum(weights * (nums - average) ** 2) / np.sum(weights))
        stds.append(std)
        averages.append(average)

    means = np.linspace(0, 18, 19 * 5)
    f = []
    for i_mean in means:
        probs = []
        for n in nums:
            prob = np.exp(-1 * i_mean) * (i_mean ** n) / np.math.factorial(n) * averages[n]
            probs.append(prob)
        f.append(np.sum(probs))

    f = np.array(f)
    my_list = np.where(means < 5)[0]

    t_popt, t_pcov = curve_fit(f=quadratic, xdata=f[my_list], ydata=means[my_list])
    t_perr = np.sqrt(np.diag(t_pcov))
    return t_popt, t_pcov, t_perr


def root_fit(x, y, yerr, model):
    # guess = array('d', [1.52157116e-11, 6.84547311e-02, 1.13069872e+07])
    guess = array('d', [(y[-1] - y[-50]) / ((x[-1] - x[-50]) * 3600 * 24) / 10000, y[0], 300 * 3600 * 24])
    names = ["p0", "p1", "L"]
    graph = ROOT.TGraphErrors(len(y), array('d', [i for i in x]), array('d', [i for i in y]),
                              array('d', [1 for i in range(len(x))]), array('d', [i for i in yerr]))

    fit = ROOT.TF1("func", model, float(x[0]), float(x[-1]), len(guess))
    for i in range(len(guess)):
        fit.SetParameter(i, guess[i])
        fit.SetParName(i, names[i])
        fit.SetParLimits(i, guess[i] - abs(guess[i]) / 2, guess[i] + abs(guess[i]) / 2)

    r = ROOT.TFitResultPtr(graph.Fit("func", "QS"))
    r.Print()
    pars = []
    errs = []
    chi = fit.GetChisquare() / fit.GetNDF()
    ndof = fit.GetNDF()

    for i in range(len(guess)):
        pars.append(fit.GetParameter(i))
        errs.append(fit.GetParError(i))

    # print(len(pars))
    # t_pcov = ROOT.TMatrixDSym(r.GetCovarianceMatrix()).GetMatrixArray()
    t_pcov = r.GetCovarianceMatrix().GetMatrixArray()
    # fitter = ROOT.TVirtualFitter.GetFitter(ROOT.TVirtualFitter())
    # print(len(t_pcov))
    pcov = np.array([0 for i in range(int(len(pars)**2))]).reshape(len(pars), len(pars))
    # print(pcov)

    '''for row in range(len(pars)):
        for col in range(len(pars)):
            pcov[row][col] = t_pcov[row][col]'''

    return pars, errs, chi, ndof


def extrapolate(model, pars, limit, typ):
    x_i = 0
    y_i = 0
    while y_i <= limit:
        y_i = model.func(np.array([x_i]), pars)[0]
        x_i += 1
    print(">", model.name, typ, limit, '{} days'.format(x_i))


def get_baseline(waveform, pre_trigger):
    baseline = 0
    for i in range(pre_trigger):
        baseline += waveform[i]
    baseline = baseline/pre_trigger
    return baseline


def get_amplitude(waveform, baseline):
    amp = 100000
    for i in range(len(waveform)):
        if waveform[i] < amp:
            amp = waveform[i]
    return amp - baseline


def get_fwhm_timestamp(waveform, baseline, peak, amplitude):
    fwhm = 0
    start = 0
    end = 0
    for i in range(peak - 50, peak + 100):
        if start == end == 0:
            if -1*(waveform[i] - baseline) > amplitude/2:
                start = i
        elif start != 0 and end == 0:
            if -1*(waveform[i] - baseline) < amplitude/2:
                end = i
        else:
            break
    fwhm = end - start
    return fwhm


def get_peak(waveform):
    peak = 0
    amp = 10000
    for i in range(len(waveform)):
        if waveform[i] < amp:
            amp = waveform[i]
            peak = i
    return peak


def mf_waveform(waveform, template):
    norm_temp = np.sqrt(np.dot(template, template))
    norm_wav = np.sqrt(np.dot(waveform, waveform))
    amplitude = np.dot(template, waveform)/norm_temp
    shape = amplitude/norm_wav

    return shape, amplitude


def lorentzian(x, amp1, cen1, wid1):
    return amp1 * wid1 ** 2 / ((x - cen1) ** 2 + wid1 ** 2)

