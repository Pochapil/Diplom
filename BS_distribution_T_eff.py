import numpy as np
from scipy.integrate import odeint
import scipy.special as special

import config  # const


def get_Teff_distribution(R_e, delta_ns, A_normal):
    # решение зависит от n размера пространства !!! взял n=3 везде;

    # над 10 формулой:
    # The case n = 2 corresponds to the spherically diverging accretion column;
    # the case n = 3 is a good approximation to the flow pattern near the magnetic poles
    # when matter falls along the field lines of a magnetic dipole.

    l0 = A_normal / (2 * delta_ns)  # длина аккреции на поверхности, м. взял как в статье
    u0 = 3 * config.H ** 2 / 8 / np.pi  # значение плотности излучения на поверхности
    d0 = delta_ns  # ширина аккреции на поверхности м
    # the amount of matter s falling on to a unit surface area of the neutron star per unit time as a fixed quantity:
    s = config.M_accretion_rate / A_normal

    # 50 формула статья
    gamma = (config.c * config.R_ns * A_normal * 3) / \
            (config.k * delta_ns ** 2 * config.M_accretion_rate * 2 * config.ksi_rad)
    # print("gamma = %f" % gamma)

    # 51 формула статья
    eta = ((8 * config.k * u0 * delta_ns ** 2 * 2 * config.ksi_rad) /
           (21 * config.c * (2 * config.G * config.M_ns * config.R_ns) ** (1 / 2) * 3)) ** (1 / 4)

    # print("eta = %f" % eta)

    # 30 формула, du/dksi; dv/dksi = производная от 3 равенства
    # возвращает u, v
    def func(y, ksi, params):
        u, v = y  # unpack current values of y
        gamma, s, G, M, R = params  # unpack parameters
        derivs = [3 * s * G * M / R * ksi ** (-5) / v,  # list of dy/dt=f functions
                  gamma * v - 3 * v / ksi - 9 / 4 * s * G * M / R * ksi ** (-5) / u]
        return derivs

    # 34 формула - из нее нахожу с помощью метода ньютона
    def findKsiShock():
        def f(x):
            return eta * gamma ** (1 / 4) * x ** (7 / 8) - 1 - np.exp(gamma * x) * \
                   (x * special.expn(2, gamma) - special.expn(2, gamma * x))

        def df(x):
            return 7 / 8 * eta * gamma ** (1 / 4) * x ** (-1 / 8) - gamma * np.exp(gamma * x) * \
                   (x * special.expn(2, gamma) - special.expn(2, gamma * x)) - np.exp(gamma * x) * \
                   (special.expn(2, gamma) + gamma * special.expn(1, gamma * x))

        def nuton(x):
            return x - f(x) / df(x)

        delta = 0.001  # точность для метода ньютона
        ksi1 = 14.3  # начальное предположение
        ksi2 = nuton(ksi1)
        while np.abs((ksi1 - ksi2)) > delta:
            ksi1 = ksi2
            ksi2 = nuton(ksi1)
        return ksi2  # rs/R - находим радиус ударной волны

    ksiShock = findKsiShock()
    # 31 формула - значения функций в точке ksiShock - граничные значения для численных расчетов
    # зависит от n размера пространства !!! взял n=3 везде
    v1 = -1 / 7 * (2 * config.G * config.M_ns / config.R_ns) ** (1 / 2) * ksiShock ** (-1 / 2)
    u1 = -3 / 4 * s * (config.G * config.M_ns / config.R_ns) * ksiShock ** (-4) / v1

    # Bundle parameters for ODE solver
    params = [gamma, s, config.G, config.M_ns, config.R_ns]
    # Bundle initial conditions for ODE solver
    y0 = [u1, v1]

    ksiStop1 = 1.
    ksi_inc = (ksiShock - ksiStop1) / (config.N_theta_accretion - 1)
    ksi_range = np.array([ksiStop1 + ksi_inc * i for i in range(config.N_theta_accretion)])
    ksi1 = ksi_range[::-1]
    solution_before_ksi = odeint(func, y0, ksi1, args=(params,), mxstep=5000000)  # от 0 до ксишок

    # analytic solve bs
    # 35 формула
    # доля излучения в стороны от всей Lt полной светимости, темп аккреции
    beta = 1 - gamma * np.exp(gamma) * (special.expn(1, gamma) - special.expn(1, gamma * ksiShock))

    # 32 формула - аналитическое решение
    def u(ksi):
        return u0 * (1 - np.exp(gamma) / beta * (special.expn(2, gamma) - special.expn(2, gamma * ksi) / ksi)) ** 4

    def v(ksi):
        return (3 / 4 * s * config.G * config.M_ns / config.R_ns * np.exp(gamma * ksi) / (ksi ** 3) * (
                1 / ksi * special.expn(2, gamma * ksi) + beta * np.exp(-gamma) - special.expn(2, gamma))) / -u(ksi)

    # переворачиваю массив потому что считал с крайней точки к 1. за границей считать нет смысла - нефизично
    u_numerical_solution = solution_before_ksi[::-1, 0]

    T = (u_numerical_solution / config.a_rad_const) ** (1 / 4)
    ksi_bs = ksi1[::-1]
    Tbs = (u(ksi_bs) / config.a_rad_const) ** (1 / 4)  # настоящее аналитическое решение

    e = config.c / (config.k * s * delta_ns)  # формула 18 стр 14
    e = config.c / config.k / (config.M_accretion_rate / A_normal) / delta_ns

    # 21 стр конец 2 абзаца
    def fTheta():
        u = solution_before_ksi[::-1, 0]
        v = solution_before_ksi[::-1, 1]
        x = ksi1[::-1]
        return -2 / 3 * e * x ** (3 / 2) * u * v

    # 21 стр конец 2 абзаца
    def fThetabs(x):
        return -2 / 3 * e * x ** (3 / 2) * u(x) * v(x)

    # ro = 1  # плотность падающего газа
    # Fr(ksi) 30 формула 17 стр
    def fr(x):
        return 4 / 3 * u(x) * v(x) + s * config.G * config.M_ns / config.R_ns * x ** (-4)

    # 19 стр под конец
    def q(ksi):
        return (ksi ** 3 * fr(ksi) - ksiShock ** 3 * fr(ksiShock)) * config.R_ns / (s * config.G * config.M_ns)

    # 30 формула 3 уравнение
    def frCalc(u, v, x):
        return 4 / 3 * u * v + s * config.G * config.M_ns / config.R_ns * x ** (-4)

    def qCalc(u, v, ksi):
        return (ksi ** 3 * frCalc(u, v, ksi) - ksiShock ** 3 * frCalc(u, v, ksi)) * config.R_ns / (
                s * config.G * config.M_ns)

    # получаем эффективную температуру из закона Стефана-Больцмана
    Teff = (fTheta() / config.sigm_Stf_Bolc) ** (1 / 4)
    Teffbs = (fThetabs(ksi_bs) / config.sigm_Stf_Bolc) ** (1 / 4)

    # формула 37, 1 - полная светимость
    L_x = (1 - beta) * config.M_accretion_rate * config.G * config.M_ns / config.R_ns

    # print("bettaBS = %f" % betta)
    # print("e = %.5f" % e)

    return Teff, ksiShock, L_x, beta
