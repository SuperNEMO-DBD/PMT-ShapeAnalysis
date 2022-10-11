"""
This file should be called in all scripts as it will then contain all functions
and routines needed for plotting etc
"""
import sys
sys.path.insert(1, '../..')

from functions.other_functions import *


"""
p_i(t) = KART/lV p_e(t - L) - 12L/pi**2  KART/lV p_e sum(1, 50, (-1)**2n/n^2 exp(-(n pi)**2 t / 6L)
p0 = KART/lV
L = L
p_i(t) = p0 p_e[(t - L - 12L/pi**2 sum(1, 50, (-1)**2n/n^2 exp(-(n pi)**2 t / 6L)] + P1

"""


class Model:
    def __init__(self):
        self.name = 'model'

    @staticmethod
    def __call__(x, pars):
        p = 101325 / 10
        p0 = pars[0] * 1E-12
        p1 = pars[1] * 1E-3
        L = pars[2] * 3600 * 24
        t = x[0] * 3600 * 24

        temp = 0
        for n in range(1, 50):
            temp += ((-1) ** n / n ** 2) * np.exp(-(n ** 2) * (np.pi ** 2) * t / (L * 6))
        f2 = (12 * L / np.pi ** 2) * temp
        y = p0 * p * (t - L - f2) + p1
        return y

    @staticmethod
    def func(x, pars):
        p = 101325 / 10
        p0 = pars[0] * 1E-12
        p1 = pars[1] * 1E-3
        L = pars[2] * 3600 * 24
        t = x * 3600 * 24

        temp = 0
        for n in range(1, 50):
            temp += ((-1) ** n / n ** 2) * np.exp(-(n ** 2) * (np.pi ** 2) * t / (L * 6))
        f2 = (12 * L / np.pi ** 2) * temp
        y = p0 * p * (t - L - f2) + p1

        return y


class Model_eff:
    def __init__(self):
        self.name = 'model_eff'

    @staticmethod
    def __call__(x, pars):
        pe = 101325 / 10
        p0 = pars[0] * 1E-12
        p1 = pars[1] * 1E-3
        L = pars[2] * 3600 * 24
        A = pars[3] * 3600 * 24
        B = pars[4] * 3600 * 24
        t_p100 = (x[0] + 98) * 3600 * 24
        t = t_p100 - (98 * 3600 * 24)

        temp = 0
        temp_p100 = 0
        for n in range(1, 50):
            temp += ((-1) ** n / n ** 2) * (np.exp(-(n ** 2) * (np.pi ** 2) * t / (L * 6)))
            temp_p100 += ((-1) ** n / n ** 2) * (np.exp(-(n ** 2) * (np.pi ** 2) * t_p100 / (L * 6)))
        f_0 = (12 * L / np.pi ** 2) * temp
        f_1 = (12 * L / np.pi ** 2) * temp_p100

        y = (p0 * pe * (t - L - f_0) + p0 * pe * 0.1 * (t_p100 - L - f_1))*(1 - 1 / (1 + np.exp(-(t - B)/A))) + p1
        return y

        '''p = 101325 / 10
        p0 = pars[0] * 1E-12
        p1 = pars[1] * 1E-3
        L = pars[2] * 3600 * 24
        A = pars[3] * 3600 * 24
        B = pars[4] * 3600 * 24
        t = x[0] * 3600 * 24

        temp = 0
        for n in range(1, 50):
            temp += ((-1) ** n / n ** 2) * (np.exp(-(n ** 2) * (np.pi ** 2) * t / (L * 6)))
        f2 = (12 * L / np.pi ** 2) * temp
        y = (p0 * p * (t - L - f2)) * (1 - 1 / (1 + np.exp(-(t - B)/A))) + p1
        return y'''

    @staticmethod
    def func(x, pars):
        pe = 101325 / 10
        p0 = pars[0] * 1E-12
        p1 = pars[1] * 1E-3
        L = pars[2] * 3600 * 24
        A = pars[3] * 3600 * 24
        B = pars[4] * 3600 * 24
        t_p100 = (x + 98) * 3600 * 24
        t = t_p100 - (98 * 3600 * 24)

        temp = 0
        temp_p100 = 0
        for n in range(1, 50):
            temp += ((-1) ** n / n ** 2) * (np.exp(-(n ** 2) * (np.pi ** 2) * t / (L * 6)))
            temp_p100 += ((-1) ** n / n ** 2) * (np.exp(-(n ** 2) * (np.pi ** 2) * t_p100 / (L * 6)))
        f_0 = (12 * L / np.pi ** 2) * temp
        f_1 = (12 * L / np.pi ** 2) * temp_p100

        y = (p0 * pe * (t - L - f_0) + p0 * pe * 0.1 * (t_p100 - L - f_1)) * (1 - 1 / (1 + np.exp(-(t - B) / A))) + p1
        return y
        '''p = 101325 / 10
        p0 = pars[0] * 1E-12
        p1 = pars[1] * 1E-3
        L = pars[2] * 3600 * 24
        A = pars[3] * 3600 * 24
        B = pars[4] * 3600 * 24
        t = x * 3600 * 24

        temp = 0
        for n in range(1, 50):
            temp += ((-1) ** n / n ** 2) * (np.exp(-(n ** 2) * (np.pi ** 2) * t / (L * 6)))
        f2 = (12 * L / np.pi ** 2) * temp
        y = (p0 * p * (t - L - f2)) * (1 - 1 / (1 + np.exp(-(t - B)/A))) + p1
        return y'''


class Model_var:
    def __init__(self):
        self.name = 'model_var'

    @staticmethod
    def __call__(x, pars, p):
        p0 = pars[0] * 1E-12
        p1 = pars[1] * 1E-3
        L = pars[2] * 3600 * 24
        t = x[0] * 3600 * 24

        temp = 0
        for n in range(1, 50):
            temp += ((-1) ** n / n ** 2) * np.exp(-(n ** 2) * (np.pi ** 2) * t / (L * 6))
        f2 = (12 * L / np.pi ** 2) * temp
        y = p0 * p * (t - L - f2) + p1

        return y

    @staticmethod
    def func(x, pars, p):
        p0 = pars[0] * 1E-12
        p1 = pars[1] * 1E-3
        L = pars[2] * 3600 * 24
        t = x * 3600 * 24

        temp = 0
        for n in range(1, 50):
            temp += ((-1) ** n / n ** 2) * np.exp(-(n ** 2) * (np.pi ** 2) * t / (L * 6))
        f2 = (12 * L / np.pi ** 2) * temp
        y = p0 * p * (t - L - f2) + p1

        return y


class Model_0:
    def __init__(self):
        self.name = 'model_0'

    @staticmethod
    def __call__(x, pars):
        pe = 101325 / 10
        p0 = pars[0] * 1E-12
        p1 = pars[1] * 1E-3
        L = pars[2] * 3600 * 24

        t_p100 = (x[0] + 98) * 3600 * 24
        t = t_p100 - (98 * 3600 * 24)
        temp = 0
        temp_p100 = 0
        for n in range(1, 50):
            temp      += ((-1) ** n / n ** 2) * (np.exp(-(n ** 2) * (np.pi ** 2) * t      / (L * 6)))
            temp_p100 += ((-1) ** n / n ** 2) * (np.exp(-(n ** 2) * (np.pi ** 2) * t_p100 / (L * 6)))
        f_0 = (12 * L / np.pi ** 2) * temp
        f_1 = (12 * L / np.pi ** 2) * temp_p100

        y = p0*pe*(t - L - f_0) + p0*pe*0.1*(t_p100 - L - f_1) + p1
        return y

    @staticmethod
    def func(x, pars):
        pe = 101325 / 10
        p0 = pars[0] * 1E-12
        p1 = pars[1] * 1E-3
        L = pars[2] * 3600 * 24

        t_p100 = (x + 98) * 3600 * 24
        t = t_p100 - (98 * 3600 * 24)
        temp = 0
        temp_p100 = 0
        for n in range(1, 50):
            temp += ((-1) ** n / n ** 2) * (np.exp(-(n ** 2) * (np.pi ** 2) * t / (L * 6)))
            temp_p100 += ((-1) ** n / n ** 2) * (np.exp(-(n ** 2) * (np.pi ** 2) * t_p100 / (L * 6)))
        f_0 = (12 * L / np.pi ** 2) * temp
        f_1 = (12 * L / np.pi ** 2) * temp_p100

        y = p0*pe*(t - L - f_0) + p0*pe*0.1*(t_p100 - L - f_1) + p1
        return y
