"""
This file should be called in all scripts as it will then contain all functions
and routines needed for plotting etc
"""
import sys
sys.path.insert(1, '../..')

from src.PMT_Classes import *
from functions.other_functions import *
from pmt_he_study.format_plot import *


class Model:
    def __init__(self):
        self.name = 'model'

    @staticmethod
    def __call__(x, pars):
        p = 101325 / 10
        p0 = pars[0]
        p1 = pars[1]
        L = pars[2]
        t = x[0] * 3600 * 24

        temp = 0
        for n in range(1, 50):
            temp += ((-1) ** n / n ** 2) * (1 - np.exp(-(n ** 2) * (np.pi ** 2) * t / (L * 6)))
        f2 = (12 / np.pi ** 2) * L * temp
        y = p0 * p * (t + f2) + p1
        return y

    @staticmethod
    def func(x, pars):
        p = 101325 / 10
        p0 = pars[0]
        p1 = pars[1]
        L = pars[2]
        t = x * 3600 * 24

        temp = 0
        for n in range(1, 50):
            temp += ((-1) ** n / n ** 2) * (1 - np.exp(-(n ** 2) * (np.pi ** 2) * t / (L * 6)))
        f2 = (12 / np.pi ** 2) * L * temp
        y = p0 * p * (t + f2) + p1
        return y


class Model_eff:
    def __init__(self):
        self.name = 'model_eff'

    @staticmethod
    def __call__(x, pars):
        p = 101325 / 10
        p0 = pars[0]
        p1 = pars[1]
        L = pars[2]
        A = pars[3] / 3600 * 24
        B = pars[4] * 3600 * 24
        t = x[0] * 3600 * 24

        temp = 0
        for n in range(1, 50):
            temp += ((-1) ** n / n ** 2) * (1 - np.exp(-(n ** 2) * (np.pi ** 2) * t / (L * 6)))
        f2 = (12 / np.pi ** 2) * L * temp
        y = (p0 * p * (t + f2) + p1) * (1 - 1 / (1 + np.exp(A * (t - B))))
        return y

    @staticmethod
    def func(x, pars):
        p = 101325 / 10
        p0 = pars[0]
        p1 = pars[1]
        L = pars[2]
        A = pars[3] / 3600 * 24
        B = pars[4] * 3600 * 24
        t = x * 3600 * 24

        temp = 0
        for n in range(1, 50):
            temp += ((-1) ** n / n ** 2) * (1 - np.exp(-(n ** 2) * (np.pi ** 2) * t / (L * 6)))
        f2 = (12 / np.pi ** 2) * L * temp
        y = p0 * p * (t + f2) + p1 * (1 - 1 / (1 + np.exp(A * (t - B))))
        return y


class Model_var:
    def __init__(self):
        self.name = 'model_var'

    @staticmethod
    def __call__(x, pars, p):
        p0 = pars[0]
        p1 = pars[1]
        L = pars[2]
        t = x[0] * 3600 * 24

        temp = 0
        for n in range(1, 50):
            temp += ((-1) ** n / n ** 2) * (1 - np.exp(-(n ** 2) * (np.pi ** 2) * t / (L * 6)))
        f2 = (12 / np.pi ** 2) * L * temp
        y = p0 * p * (t + f2) + p1
        return y

    @staticmethod
    def func(x, pars, p):
        p0 = pars[0]
        p1 = pars[1]
        L = pars[2]
        t = x * 3600 * 24

        temp = 0
        for n in range(1, 50):
            temp += ((-1) ** n / n ** 2) * (1 - np.exp(-(n ** 2) * (np.pi ** 2) * t / (L * 6)))
        f2 = (12 / np.pi ** 2) * L * temp
        y = p0 * p * (t + f2) + p1
        return y


class Model_0:
    def __init__(self):
        self.name = 'model_0'

    @staticmethod
    def __call__(x, pars):
        pe_0 = 101325 / 100
        pe_1 = 101325 / 10
        p0 = pars[0]
        p1 = pars[1]
        L = pars[2]

        t_0 = (x[0] + 98) * 3600 * 24
        t_1 = t_0 - (98 * 3600 * 24)
        temp_0 = 0
        temp_1 = 0
        for n in range(1, 50):
            temp_0 += ((-1) ** n / n ** 2) * (1 - np.exp(-(n ** 2) * (np.pi ** 2) * t_0 / (L * 6)))
            temp_1 += ((-1) ** n / n ** 2) * (1 - np.exp(-(n ** 2) * (np.pi ** 2) * t_1 / (L * 6)))
        f_0 = (12 / np.pi ** 2) * L * temp_0
        f_1 = (12 / np.pi ** 2) * L * temp_1

        y = p0 * (pe_0 * (t_0 + f_0) + pe_1 * (t_1 + f_1)) + p1
        return y
    @staticmethod
    def func(x, pars):
        pe_0 = 101325 / 100
        pe_1 = 101325 / 10
        p0 = pars[0]
        p1 = pars[1]
        L = pars[2]

        t_0 = (x + 98) * 3600 * 24
        t_1 = t_0 - (98 * 3600 * 24)
        temp_0 = 0
        temp_1 = 0
        for n in range(1, 50):
            temp_0 += ((-1) ** n / n ** 2) * (1 - np.exp(-(n ** 2) * (np.pi ** 2) * t_0 / (L * 6)))
            temp_1 += ((-1) ** n / n ** 2) * (1 - np.exp(-(n ** 2) * (np.pi ** 2) * t_1 / (L * 6)))
        f_0 = (12 / np.pi ** 2) * L * temp_0
        f_1 = (12 / np.pi ** 2) * L * temp_1

        y = p0 * (pe_0 * (t_0 + f_0) + pe_1 * (t_1 + f_1)) + p1
        return y
