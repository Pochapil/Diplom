import numpy as np
from scipy.integrate import odeint
import scipy.special as special
import matplotlib.pyplot as plt

import config  # const


# 34 формула - из нее нахожу с помощью метода ньютона
def find_ksi_shock(eta, gamma):
    '''Методом Ньютона нахожу'''

    def f(x):
        return eta * gamma ** (1 / 4) * x ** (7 / 8) - 1 - np.exp(gamma * x) * \
               (x * special.expn(2, gamma) - special.expn(2, gamma * x))

    def df(x):
        # df/dx
        return 7 / 8 * eta * gamma ** (1 / 4) * x ** (-1 / 8) - gamma * np.exp(gamma * x) * \
               (x * special.expn(2, gamma) - special.expn(2, gamma * x)) - np.exp(gamma * x) * \
               (special.expn(2, gamma) + gamma * special.expn(1, gamma * x))

    def nuton(x):
        return x - f(x) / df(x)

    delta = 0.001  # точность для метода ньютона
    ksi_prev = 31  # начальное предположение
    ksi_next = nuton(ksi_prev)
    while np.abs((ksi_prev - ksi_next)) > delta:
        ksi_prev = ksi_next
        ksi_next = nuton(ksi_prev)
    return ksi_next  # rs/R - находим радиус ударной волны


def solve_numerical(gamma, s, ksi_shock):
    '''Численно решить и найти распределение u,v '''
    # 30 формула, du/dksi; dv/dksi = производная от 3 равенства
    # возвращает u, v
    # dv/dksi du/dksi
    '''n = 3. подставил Fr во 2 выражение. раскрыл производную и выразил dv/dksi. 
    du/dksi выражается через 1 уравнение. Это функции для производных'''

    def func(y, ksi, params):
        u, v = y  # unpack current values of y
        gamma, s, G, M, R = params  # unpack parameters
        derivs = [3 * s * G * M / R * ksi ** (-5) / v,  # list of dy/dt=f functions
                  gamma * v - 3 * v / ksi - 9 / 4 * s * G * M / R * ksi ** (-5) / u]
        return derivs

    # 31 формула - значения функций в точке ksiShock - граничные значения для численных расчетов
    # зависит от n размера пространства !!! взял n=3 везде
    '''граничные условия'''
    v1 = -1 / 7 * (2 * config.G * config.M_ns / config.R_ns) ** (1 / 2) * ksi_shock ** (-1 / 2)
    u1 = -3 / 4 * s * (config.G * config.M_ns / config.R_ns) * ksi_shock ** (-4) / v1

    # Bundle parameters for ODE solver
    params = [gamma, s, config.G, config.M_ns, config.R_ns]
    # Bundle initial conditions for ODE solver
    y0 = [u1, v1]

    ksi_start = 1.
    ksi_range = np.linspace(ksi_start, ksi_shock, config.N_theta_accretion)
    # переворачиваем так как идем от ударной волны ???
    ksi_range = ksi_range[::-1]
    solution_before_ksi = odeint(func, y0, ksi_range, args=(params,), mxstep=5000000)  # от 0 до ксишок

    # переворачиваю массив потому что считал с крайней точки к 1. за границей считать нет смысла - нефизично
    u_numerical_solution = solution_before_ksi[::-1, 0]
    v_numerical_solution = solution_before_ksi[::-1, 1]
    return u_numerical_solution, v_numerical_solution


# 32 формула - аналитическое решение
def get_u_analytic(u0, gamma, beta, ksi):
    return u0 * (1 - np.exp(gamma) / beta * (special.expn(2, gamma) - special.expn(2, gamma * ksi) / ksi)) ** 4


def get_v_analytic(u0, s, gamma, beta, ksi):
    u = get_u_analytic(u0, gamma, beta, ksi)
    return (3 / 4 * s * config.G * config.M_ns / config.R_ns * np.exp(gamma * ksi) / (ksi ** 3) * (
            1 / ksi * special.expn(2, gamma * ksi) + beta * np.exp(-gamma) - special.expn(2, gamma))) / -u


# 21 стр конец 2 абзаца
def get_f_theta(u, v, ksi_arr, e):
    # Numerical
    # f_theta is the energy flux radiated by the unit area of this surface
    return -2 / 3 * e * ksi_arr ** (3 / 2) * u * v


# 21 стр конец 2 абзаца
def get_f_theta_bs(ksi_arr, e, u0, s, gamma, beta):
    # analytical
    return -2 / 3 * e * ksi_arr ** (3 / 2) * get_u_analytic(u0, gamma, beta, ksi_arr) * get_v_analytic(u0, s, gamma,
                                                                                                       beta, ksi_arr)


# ro = 1 # плотность падающего газа
# Fr(ksi) 30 формула 17 стр
def fr(s, u0, gamma, beta, ksi_arr):
    return 4 / 3 * get_u_analytic(u0, gamma, beta, ksi_arr) * get_v_analytic(u0, s, gamma, beta, ksi_arr) + \
           s * config.G * config.M_ns / config.R_ns * ksi_arr ** (-4)


# 19 стр под конец
def q(ksi, ksi_shock, s):
    return (ksi ** 3 * fr(ksi) - ksi_shock ** 3 * fr(ksi_shock)) * config.R_ns / (s * config.G * config.M_ns)


# 30 формула 3 уравнение
def frCalc(u, v, x, s):
    return 4 / 3 * u * v + s * config.G * config.M_ns / config.R_ns * x ** (-4)


def qCalc(u, v, ksi, s, ksi_shock):
    return (ksi ** 3 * frCalc(u, v, ksi, s) - ksi_shock ** 3 * frCalc(u, v, ksi, s)) * config.R_ns / (
            s * config.G * config.M_ns)


def get_Teff_distribution(delta_ns, A_normal):
    '''решение зависит от n размера пространства !!! взял n=3 везде
    delta_ns - уповерхности НЗ --- the distance between the field lines $\delta$
    A_normal - уповерхности НЗ --- a cross-section
    A_perp =  2, delta, 2 uppi a R sin theta --- 2!!! это площадка для 2 колонок!!
    '''

    # над 10 формулой:
    # The case n = 2 corresponds to the spherically diverging accretion column;
    # the case n = 3 is a good approximation to the flow pattern near the magnetic poles
    # when matter falls along the field lines of a magnetic dipole.

    # l0 = A_normal / (2 * delta_ns)  # длина аккреции на поверхности, м. взял как в статье
    # d0 = delta_ns  # ширина аккреции на поверхности м

    u0 = 3 * config.H ** 2 / 8 / np.pi  # значение плотности излучения на поверхности
    # the amount of matter s falling on to a unit surface area of the neutron star per unit time as a fixed quantity:
    s = config.M_accretion_rate / A_normal

    # 50 формула статья Dimensionless coefficients
    gamma = (config.c * config.R_ns * A_normal * 3) / \
            (config.k * delta_ns ** 2 * config.M_accretion_rate * 2 * config.ksi_rad)

    # 51 формула статья
    eta = ((8 * config.k * u0 * delta_ns ** 2 * 2 * config.ksi_rad) /
           (21 * config.c * (2 * config.G * config.M_ns * config.R_ns) ** (1 / 2) * 3)) ** (1 / 4)

    e = config.c / (config.k * s * delta_ns)  # формула 18 стр 14

    ksi_shock = find_ksi_shock(eta, gamma)
    # print(ksiShock)

    ksi_arr = np.linspace(1., ksi_shock, config.N_theta_accretion)

    numerical_flag = False
    if numerical_flag:
        u_numerical_solution, v_numerical_solution = solve_numerical(gamma, s, ksi_shock)
        T = (u_numerical_solution / config.a_rad_const) ** (1 / 4)
        Teff = (get_f_theta(u_numerical_solution, v_numerical_solution, ksi_arr, e) / config.sigm_Stf_Bolc) ** (1 / 4)

    # 35 формула
    # доля излучения в стороны от всей Lt полной светимости, темп аккреции
    beta = 1 - gamma * np.exp(gamma) * (special.expn(1, gamma) - special.expn(1, gamma * ksi_shock))

    # analytic solve bs
    v = get_v_analytic(u0, s, gamma, beta, ksi_arr)
    u = get_u_analytic(u0, gamma, beta, ksi_arr)

    Tbs = (get_u_analytic(u0, gamma, beta, ksi_arr) / config.a_rad_const) ** (1 / 4)  # настоящее аналитическое решение

    # получаем эффективную температуру из закона Стефана-Больцмана
    Teffbs = (get_f_theta_bs(ksi_arr, e, u0, s, gamma, beta) / config.sigm_Stf_Bolc) ** (1 / 4)

    # формула 37, 1 - полная светимость
    L_x = (1 - beta) * config.M_accretion_rate * config.G * config.M_ns / config.R_ns

    plot_flag = False
    if plot_flag:
        fig = plt.figure(figsize=(12, 6))
        ax = fig.add_subplot(111)
        ax.plot(ksi_arr, u, label='bs')
        ax.plot(ksi_arr, u_numerical_solution, label='numerical')
        ax.legend()
        plt.show()

        fig = plt.figure(figsize=(12, 6))
        ax = fig.add_subplot(111)
        ax.plot(ksi_arr, v, label='bs')
        ax.plot(ksi_arr, v_numerical_solution, label='numerical')
        ax.legend()
        plt.show()

        fig = plt.figure(figsize=(12, 6))
        ax = fig.add_subplot(111)
        ax.plot(ksi_arr, Tbs, label='bs')
        ax.plot(ksi_arr, T, label='numerical')
        ax.legend()
        plt.show()

        fig = plt.figure(figsize=(12, 6))
        ax = fig.add_subplot(111)
        ax.plot(ksi_arr, Teffbs, label='bs')
        ax.plot(ksi_arr, Teff, label='numerical')
        ax.legend()
        plt.show()

    return Teffbs, ksi_shock, L_x, beta


if __name__ == '__main__':
    delta_ns = 20434.12420796042
    A_normal = 4571886568.042245
    T_eff, ksi_shock, L_x, beta = get_Teff_distribution(delta_ns, A_normal)
