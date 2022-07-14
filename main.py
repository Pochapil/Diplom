import geometricTask.matrix as matrix
import numpy as np
import matplotlib.pyplot as plt

import scipy.integrate

import matplotlib.cm as cm
from matplotlib.colors import Normalize

import BS_distribution_T_eff as get_T_eff
import config  # const

# смотри RotationSingleColorBar.py
# const
file_count = 29
grad_to_rad = np.pi / 180

# Parameters
MSun = config.MSun  # масса молнца г
G = config.G  # гравитационная постоянная см3·с−2·г−1
M_ns = config.M_ns  # масса нз г
R_ns = config.R_ns  # радиус нз см
sigmStfBolc = config.sigmStfBolc  # постоянная стефана больцмана в сгс

# new args (for new func)
dRe_div_Re = config.dRe_div_Re  # взял просто число
M_accretion_rate = config.M_accretion_rate  # темп аккреции
print(M_accretion_rate)
ksi_rad = config.ksi_rad
H = config.H  # магнитное поле стр 19 над формулой 37
a_portion = config.a_portion  # a - в азимутальном направлении поток занимает фиксированную долю a полного круга 2πR sinθ
mu = config.mu  # магнитный момент Гаусс * см3
ksi_param = 0.5  # между 1 и 2 формулой в статье
R_alfven = (mu ** 2 / (2 * M_accretion_rate * (2 * G * M_ns) ** (1 / 2))) ** (2 / 7)  # альфвеновский радиус
R_e = ksi_param * R_alfven  # между 1 и 2 формулой в статье
print("R_e :%f" % R_e)
print("R_e/R_* = %f" % (R_e / R_ns))

# количество шагов
N_phi_accretion = 100
N_theta_accretion = 100

# угол между нормалью к двойной системе и наблюдателем
i_angle = 0 * grad_to_rad
# вектор на наблюдателя в системе координат двойной системы
e_obs = np.array([0, np.sin(i_angle), np.cos(i_angle)])

# угол между осью вращения системы и собственным вращенеим НЗ
betta_rotate = (file_count // 3) * 15 * grad_to_rad
betta_rotate = 60 * grad_to_rad
phi_rotate = 0 * grad_to_rad

# угол между собственным вращенеим НЗ и магнитной осью
betta_mu = (file_count % 3) * 15 * grad_to_rad
betta_mu = 80 * grad_to_rad
phi_mu_0 = 0 * grad_to_rad


# формула 2 в статье
def delta_distance(theta):
    # R=R_e * sin_theta ** 2
    return R_e * np.sin(theta) ** 3 / (1 + 3 * np.cos(theta) ** 2) ** (1 / 2) * dRe_div_Re


# формула 3 в статье
def A_normal(theta):
    # a - в азимутальном направлении поток занимает фиксированную долю a полного круга 2πR sinθ
    return 2 * delta_distance(theta) * 2 * np.pi * a_portion * R_e * np.sin(theta) ** 3


theta_accretion_begin = np.arcsin((R_ns / R_e) ** (1 / 2))  # от поверхности NS - угол при котором радиус = радиусу НЗ

print("A|(R*)/R*2 = %f" % (A_normal(theta_accretion_begin) / R_ns ** 2))
print("delta(R*)/R* = %f" % (delta_distance(theta_accretion_begin) / R_ns))

Teff, ksiShock, L_x = get_T_eff.get_Teff_distribution(N_theta_accretion, R_e, delta_distance(theta_accretion_begin),
                                                      A_normal(theta_accretion_begin))
# ksiShock = 4.523317
print("ksiShock: %f" % ksiShock)
print("Rshock/R*: %f" % ksiShock)
print("table 2 column 3  %f" % (G * M_ns * M_accretion_rate / (R_ns * config.L_edd)))  # значение в табличке
# print("arcsin begin: %f" % (R_ns / R_e) ** (1 / 2))
# print("arcsin end: %f" % (R_ns * ksiShock / R_e) ** (1 / 2))

theta_accretion_end = np.arcsin((R_ns * ksiShock / R_e) ** (1 / 2))  # до шока - угол когда радиус = радиус шока

R_phi = R_ns * np.sin(theta_accretion_begin)

l0 = A_normal(theta_accretion_begin) / (2 * delta_distance(theta_accretion_begin))
print("l0 = %.5f" % l0)
lim_phi_accretion = 360 * grad_to_rad * 1.01 * a_portion  # полный круг для наложения. на карте были пробелы

print("phi = %f" % (lim_phi_accretion / grad_to_rad))
print("theta_accretion_begin = %f" % (theta_accretion_begin / grad_to_rad))
print("theta_accretion_end = %f" % (theta_accretion_end / grad_to_rad))

# шаги по углам для интегрирования
step_phi_accretion = lim_phi_accretion / N_phi_accretion
step_theta_accretion = (theta_accretion_end - theta_accretion_begin) / N_theta_accretion
print("step theta = %f " % step_theta_accretion)

# для отрисовки карты и интеграла. создаю пустые массивы для распределения углов и косинуса
phi_range = np.array([step_phi_accretion * i for i in range(N_phi_accretion)])
theta_range = np.array([theta_accretion_begin + step_theta_accretion * j for j in range(N_theta_accretion)])
cos_psi_range = np.empty([N_phi_accretion, N_theta_accretion])

# тут беру значения из Teff в промежутках так как нахожу по отрезку кси а нужно по тета!!

ksiStop1 = 1.
ksiInc1 = - (ksiShock - ksiStop1) / N_theta_accretion
ksi1 = np.arange(ksiShock, ksiStop1, ksiInc1)
ksi_bs = ksi1[::-1]
i = 0  # left_border
true_T_eff = []
for theta in theta_range:
    while ksi_bs[i + 1] < R_e / R_ns * np.sin(theta) ** 2 and (i < len(theta_range) - 2):
        i += 1
    true_T_eff.append(Teff[i])

fig = plt.figure(figsize=(8, 8))
ax1 = fig.add_subplot(111)
ax1.plot(ksi_bs, Teff, 'b')
ax1.plot(ksi_bs, true_T_eff, 'r')
Teff = np.array(true_T_eff)


def create_array_normal(phi_range, theta_range):
    N_phi_accretion = len(phi_range)
    N_theta_accretion = len(theta_range)
    array_normal = []  # матрица нормалей чтобы не пересчитывать в циклах
    # matrix_normal = np.empty([N_phi_accretion, N_theta_accretion])
    for i in range(N_phi_accretion):
        for j in range(N_theta_accretion):
            # matrix_normal[i, j] = matrix.newE_n(phi_range[i], theta_range[j])
            array_normal.append(matrix.newE_n(phi_range[i], theta_range[j]))
    return array_normal


array_normal = create_array_normal(phi_range, theta_range)

dS = []  # массив единичных площадок при интегрировании так как зависит только от theta посчитаю 1 раз
dS_simps = []
# формула 5 из статьи для dl
for j in range(N_theta_accretion):
    dl = R_e * (3 * np.cos(theta_range[j]) ** 2 + 1) ** (1 / 2) * np.sin(theta_range[j]) * step_theta_accretion
    dphi = R_e * np.sin(theta_range[j]) ** 3 * step_phi_accretion  # R=R_e * sin_theta ** 2; R_phi = R * sin_theta
    dS.append(dphi * dl)  # единичная площадка при интегрировании
    # аналогично для интеграла по симпсону
    dl_simps = R_e * (3 * np.cos(theta_range[j]) ** 2 + 1) ** (1 / 2) * np.sin(theta_range[j])
    dphi_simps = R_e * np.sin(theta_range[j]) ** 3
    dS_simps.append(dphi_simps * dl_simps)

omega_ns = 8  # скорость вращения НЗ - будет меняться только угол phi_mu!
# цикл для поворотов, сколько точек на графике интегралов - для фазы от 0 до 2 - с перекрытием чтобы форму макс
max_phase_angle = 720
t_max = (max_phase_angle // omega_ns) + (1 if max_phase_angle % omega_ns > 0 else 0)
# цикл для поворотов, сколько точек на графике интегралов - для фазы от 0 до 1 (полного поворота)
max_phase_angle = 360
t_max_for_cos = (max_phase_angle // omega_ns) + (1 if max_phase_angle % omega_ns > 0 else 0)

omega_ns = omega_ns * grad_to_rad  # перевожу в радианы


def get_angles_from_vector(vector):
    x = vector[0, 0]
    y = vector[0, 1]
    z = vector[0, 2]
    theta = np.arccos(z)  # np.arccos(z/r)
    if x > 0:
        if y >= 0:
            phi = np.arctan(y / x)
        else:
            phi = np.arctan(y / x) + 2 * np.pi
    else:
        phi = np.arctan(y / x) + np.pi
    return phi, theta


def get_angles_from_vector_one_dimension(vector):
    x = vector[0]
    y = vector[1]
    z = vector[2]
    theta = np.arccos(z)  # np.arccos(z/r)
    if x > 0:
        if y >= 0:
            phi = np.arctan(y / x)
        else:
            phi = np.arctan(y / x) + 2 * np.pi
    else:
        phi = np.arctan(y / x) + np.pi
    return phi, theta


def get_lim_for_analytic_integral_phi(theta, e_obs):
    # мб нужно поднять чтобы не считать углы наблюдателя, а посчитать 1 раз
    phi_obs, theta_obs = get_angles_from_vector(e_obs)

    # вроде в тетрадке есть решение
    def get_limit_delta_phi(theta, theta_obs):
        divider = (1 - 3 * np.cos(theta) ** 2) * np.sin(theta_obs)
        lim = - 3 * np.sin(theta) * np.cos(theta) * np.cos(theta_obs) / divider
        if divider > 0:
            lim = -lim
        # (phi_range[i] >= delta_phi_lim) and (phi_range[i] <= 2 * np.pi - delta_phi_lim) - условие
        if lim >= 1:
            return 0  # любой угол будет больше 0 и меньше 2 * np.pi
        if lim <= -1:
            # чтобы интеграл дал 0 нужно pi: от pi до 2 pi - pi
            return np.pi  # все углы будут меньше 2 * np.pi и 1 скобка даст false
        return np.arccos(lim)

    return get_limit_delta_phi(theta, theta_obs)  # arccos < delta_phi < 2 pi - arccos


def calculate_total_luminosity(phi_range, theta_range):
    N_phi_accretion = len(phi_range)
    total_luminosity_step = [0] * N_phi_accretion
    for i in range(N_phi_accretion):
        total_luminosity_step[i] = sigmStfBolc * scipy.integrate.simps(Teff ** 4 * dS_simps, theta_range)
    total_luminosity = scipy.integrate.simps(total_luminosity_step, phi_range)
    return total_luminosity


def check_if_intersect(origin_phi, origin_theta, direction_vector, lim_phi_accretion, lim_theta):
    # все отнормирую на радиус НЗ
    # r = origin + t * direction - уравнение луча

    x_origin = np.sin(origin_theta) * np.cos(origin_phi)
    y_origin = np.sin(origin_theta) * np.sin(origin_phi)
    z_origin = np.cos(origin_theta)
    # origin_point = [x, y, z]

    x_direction = direction_vector[0, 0]
    y_direction = direction_vector[0, 1]
    z_direction = direction_vector[0, 2]

    def find_intersect_solution(a, b, c):
        t_1 = (-b + (b ** 2 - 4 * a * c) ** (1 / 2)) / (2 * a)
        t_2 = (-b - (b ** 2 - 4 * a * c) ** (1 / 2)) / (2 * a)
        return t_1, t_2

    # sphere x**2 + y**2 + z**2 == 1
    a_sphere = x_direction ** 2 + y_direction ** 2 + z_direction ** 2
    b_sphere = 2 * (x_origin * x_direction + y_origin * y_direction + z_origin * z_direction)
    c_sphere = x_origin ** 2 + y_origin ** 2 + z_origin ** 2 - 1

    # solution for sphere
    t_sphere = find_intersect_solution(a_sphere, b_sphere, c_sphere)

    # print("t_sphere1 = %f,t_sphere2 = %f" % (t_sphere[0], t_sphere[1]))
    if t_sphere[0] > 0 or t_sphere[1] > 0:
        return True

    # cone x**2 + y**2 == z**2
    a_cone = x_direction ** 2 + y_direction ** 2 - z_direction ** 2
    b_cone = 2 * (x_origin * x_direction + y_origin * y_direction - z_origin * z_direction)
    c_cone = x_origin ** 2 + y_origin ** 2 - z_origin ** 2

    # solution for cone
    t_cone = find_intersect_solution(a_cone, b_cone, c_cone)
    # print("t_cone1 = %f,t_cone2 = %f" % (t_cone[0], t_cone[1]))

    # для верхнего конуса:
    if t_cone[0] > 0:
        intersect_vector = np.array([x_origin, y_origin, z_origin]) + t_cone[0] * np.array(
            [x_direction, y_direction, z_direction])
        phi_intersect, theta_intersect = get_angles_from_vector_one_dimension(intersect_vector)
        if (intersect_vector[2] > 0 and theta_intersect < lim_theta and phi_intersect < lim_phi_accretion):
            return True

    if t_cone[1] > 0:
        intersect_vector = np.array([x_origin, y_origin, z_origin]) + t_cone[1] * np.array(
            [x_direction, y_direction, z_direction])
        phi_intersect, theta_intersect = get_angles_from_vector_one_dimension(intersect_vector)
        if (intersect_vector[2] > 0 and theta_intersect < lim_theta and phi_intersect < lim_phi_accretion):
            return True

    return False


def calculate_integral_distribution(t_max, N_phi_accretion, N_theta_accretion):
    integral_max = 0
    # sum_intense изотропная светимость ( * 4 pi еще надо)
    sum_intense = [0] * t_max

    # для интеграла по simpson
    sum_simps_integrate = [0] * t_max
    simps_integrate_step = [0] * N_phi_accretion
    simps_cos = [0] * N_theta_accretion  # cos для интеграла по симпсону

    # для аналитического интеграла
    lim_phi_begin = [0] * N_theta_accretion
    analytic_integral_phi = [0] * t_max
    for i1 in range(t_max):
        # поворот
        phi_mu = phi_mu_0 + omega_ns * i1
        # расчет матрицы поворота в магнитную СК и вектора на наблюдателя
        A_matrix_analytic = matrix.newMatrixAnalytic(phi_rotate, betta_rotate, phi_mu, betta_mu)

        # print("analytic matrix:")
        # print(A_matrix_analytic)
        e_obs_mu = np.dot(A_matrix_analytic, e_obs)  # переход в магнитную СК
        # print("e_obs_mu%d: (%f, %f, %f), angle phi = %f" % (
        # i1, e_obs_mu[0, 0], e_obs_mu[0, 1], e_obs_mu[0, 2], np.arctan(e_obs_mu[0, 1] / e_obs_mu[0, 0])/grad_to_rad))
        # print("e_obs_mu%d: (%f, %f, %f)" % (i1, np.take(e_obs_mu, 0), np.take(e_obs_mu, 1), np.take(e_obs_mu, 2)))

        phi, theta = get_angles_from_vector(e_obs_mu)
        # print("thetaObs%d = %f" % (i1, (theta / grad_to_rad)))

        # sum_intense изотропная светимость ( * 4 pi еще надо)
        for i in range(N_phi_accretion):
            for j in range(N_theta_accretion):
                # cos_psi_range[i][j] = np.dot(e_obs_mu, matrix.newE_n(phi_range[i], theta_range[j])) # неэффективно
                cos_psi_range[i, j] = np.dot(e_obs_mu, array_normal[i * N_theta_accretion + j])

                if check_if_intersect(phi_range[i], theta_range[j], e_obs_mu, lim_phi_accretion, theta_accretion_end):
                    simps_cos[j] = 0

                elif cos_psi_range[i, j] > 0:
                    sum_intense[i1] += sigmStfBolc * Teff[j] ** 4 * cos_psi_range[i, j] * dS[j]
                    simps_cos[j] = cos_psi_range[i, j]
                    # * S=R**2 * step_phi_accretion * step_theta_accretion
                else:
                    simps_cos[j] = 0

            simps_integrate_step[i] = sigmStfBolc * scipy.integrate.simps(Teff ** 4 * simps_cos * dS_simps, theta_range)
        # находим позицию максимума
        if integral_max < sum_intense[i1]:
            position_of_max = i1
            integral_max = sum_intense[i1]

        sum_simps_integrate[i1] = scipy.integrate.simps(simps_integrate_step, phi_range)

        for j in range(N_theta_accretion):
            lim_phi_begin[j] = get_lim_for_analytic_integral_phi(theta_range[j],
                                                                 e_obs_mu)  # считаем границы для интеграла

        phi_obs, theta_obs = get_angles_from_vector(e_obs_mu)

        L1 = (1 - 3 * np.array(np.cos(theta_range) ** 2)) * np.array(np.sin(theta_obs)) * (
                np.array(np.sin(2 * np.pi - np.array(lim_phi_begin))) - np.array(np.sin(lim_phi_begin))) + 3 * np.sin(
            theta_range) * np.cos(theta_range) * np.cos(theta_obs) * 2 * (np.pi - np.array(lim_phi_begin))

        L = sigmStfBolc * Teff ** 4 * R_e ** 2 * np.sin(theta_range) ** 4 * L1
        analytic_integral_phi[i1] = scipy.integrate.simps(L, theta_range)

    return sum_intense, sum_simps_integrate, analytic_integral_phi, position_of_max


sum_intense, sum_simps_integrate, analytic_integral_phi, position_of_max = calculate_integral_distribution(t_max,
                                                                                                           N_phi_accretion,
                                                                                                           N_theta_accretion)

print("max: %d" % position_of_max)


def plot_map_cos(n_pos, position_of_max, t_max, N_phi_accretion, N_theta_accretion, row_number, column_number):
    number_of_plots = row_number * column_number

    crf = [0] * number_of_plots
    cr = [0] * number_of_plots

    fig, axes = plt.subplots(nrows=row_number, ncols=column_number, figsize=(8, 8), subplot_kw={'projection': 'polar'})
    # сдвигаем графики относительно позиции максимума. чтобы макс был на (0,0)
    row_figure = 0
    column_figure = 0
    # phi_mu = phi_mu_0
    for i1 in range(number_of_plots):
        phi_mu = phi_mu_0 + omega_ns * (n_pos + i1 + position_of_max)
        # расчет матрицы поворота в магнитную СК и вектора на наблюдателя
        A_matrix_analytic = matrix.newMatrixAnalytic(phi_rotate, betta_rotate, phi_mu, betta_mu)

        # A_matrix_analytic = matrix.newRy(betta_mu) @ matrix.newRz(phi_mu) @ matrix.newRy(betta_rotate) \
        #                 @ matrix.newRz(phi_rotate)

        # print("analytic matrix:")
        # print(A_matrix_analytic)
        count_0 = 0
        e_obs_mu = np.dot(A_matrix_analytic, e_obs)  # переход в магнитную СК
        for i in range(N_phi_accretion):
            for j in range(N_theta_accretion):
                cos_psi_range[i][j] = np.dot(e_obs_mu, array_normal[i * N_phi_accretion + j])
                if cos_psi_range[i][j] < 0:
                    count_0 += 1

        crf[i1] = axes[row_figure, column_figure].contourf(phi_range, theta_range / grad_to_rad,
                                                           cos_psi_range.transpose(), vmin=-1, vmax=1)

        if count_0 > 0:
            cr[i1] = axes[row_figure, column_figure].contour(phi_range, theta_range / grad_to_rad,
                                                             cos_psi_range.transpose(),
                                                             [0.], colors='w')
        column_figure += 1
        if column_figure == column_number:
            column_figure = 0
            row_figure += 1

    cbar = fig.colorbar(crf[i1], ax=axes[:], shrink=0.8, location='right')

    plt.show()


def plot_map_cos_in_range(position_of_max, t_max_for_cos, N_phi_accretion, N_theta_accretion, row_number,
                          column_number):
    number_of_plots = row_number * column_number

    crf = [0] * number_of_plots
    cr = [0] * number_of_plots

    fig, axes = plt.subplots(nrows=row_number, ncols=column_number, figsize=(8, 8), subplot_kw={'projection': 'polar'})

    row_figure = 0
    column_figure = 0

    # для единого колорбара - взял из инета
    cmap = cm.get_cmap('viridis')
    normalizer = Normalize(-1, 1)
    im = cm.ScalarMappable(norm=normalizer)
    # сдвигаем графики относительно позиции максимума. чтобы макс был на (0,0)
    phi_mu_max = phi_mu_0 + omega_ns * position_of_max
    for i1 in range(number_of_plots):
        # поворот на угол
        phi_mu = phi_mu_max + omega_ns * (t_max_for_cos / (number_of_plots - 1)) * i1
        # расчет матрицы поворота в магнитную СК и вектора на наблюдателя
        A_matrix_analytic = matrix.newMatrixAnalytic(phi_rotate, betta_rotate, phi_mu, betta_mu)

        count_0 = 0  # счетчик для контура 0 на карте
        e_obs_mu = np.dot(A_matrix_analytic, e_obs)  # переход в магнитную СК
        # e_obs_mu = np.array([0,1,-1])
        for i in range(N_phi_accretion):
            for j in range(N_theta_accretion):
                cos_psi_range[i][j] = np.dot(e_obs_mu, array_normal[i * N_phi_accretion + j])
                if cos_psi_range[i][j] < 0:
                    count_0 += 1

        crf[i1] = axes[row_figure, column_figure].contourf(phi_range, theta_range / grad_to_rad,
                                                           cos_psi_range.transpose(), vmin=-1, vmax=1, cmap=cmap,
                                                           norm=normalizer)
        # попытки для сдвига 0 на картах
        # axes[row_figure, column_figure].set_ylim(theta_range[0]/grad_to_rad, theta_range[N_theta_accretion-1]/grad_to_rad)
        axes[row_figure, column_figure].set_rorigin(-theta_accretion_begin)  # отступ для центра звезды

        axes[row_figure, column_figure].set_yticks([(theta_range[0] / grad_to_rad),
                                                    (theta_range[N_theta_accretion // 2] / grad_to_rad),
                                                    (theta_range[-1] / grad_to_rad)])  # форматируем по тета

        axes[row_figure, column_figure].set_yticklabels([round(theta_range[0] / grad_to_rad, 1),
                                                         round(theta_range[N_theta_accretion // 2] / grad_to_rad, 1),
                                                         round(theta_range[-1] / grad_to_rad,
                                                               1)])  # форматируем по тета

        # axes[row_figure, column_figure].set_theta_zero_location('W', offset=50)
        if count_0 > 0:  # рисую контур 0 если он есть
            cr[i1] = axes[row_figure, column_figure].contour(phi_range, theta_range / grad_to_rad,
                                                             cos_psi_range.transpose(),
                                                             [0.], colors='w')

        axes[row_figure, column_figure].set_title(
            "phase = %.2f" % (omega_ns * (t_max_for_cos / (number_of_plots - 1)) * i1 / (2 * np.pi)))
        column_figure += 1
        if column_figure == column_number:
            column_figure = 0
            row_figure += 1

    plt.subplots_adjust(hspace=0.5, wspace=0.5)
    cbar = fig.colorbar(im, ax=axes[:, :], shrink=0.7, location='right')
    plt.show()


def plot_map_t_eff(T_eff, N_phi_accretion, N_theta_accretion):
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    # fig = plt.figure(figsize=(8, 8), projection="polar")
    # ax = fig.add_subplot(111)
    result = np.empty([N_phi_accretion, N_theta_accretion])
    for i in range(N_phi_accretion):
        for j in range(N_theta_accretion):
            result[i, j] = T_eff[j]
    ax.contourf(phi_range, theta_range / grad_to_rad, result.transpose())
    plt.show()


def plot_map_delta_phi(position_of_max, t_max_for_cos, N_phi_accretion, N_theta_accretion, row_number, column_number):
    number_of_plots = row_number * column_number

    crf = [0] * number_of_plots

    fig, axes = plt.subplots(nrows=row_number, ncols=column_number, figsize=(8, 8), subplot_kw={'projection': 'polar'})

    row_figure = 0
    column_figure = 0

    # для единого колорбара - взял из инета
    cmap = cm.get_cmap('viridis')
    normalizer = Normalize(-1, 1)
    im = cm.ScalarMappable(norm=normalizer)

    phi_mu_max = phi_mu_0 + omega_ns * position_of_max
    for i1 in range(number_of_plots):
        # сдвигаем графики относительно позиции максимума. чтобы макс был на (0,0)
        phi_mu = phi_mu_max + omega_ns * (t_max_for_cos / (number_of_plots - 1)) * i1
        # расчет матрицы поворота в магнитную СК и вектора на наблюдателя
        A_matrix_analytic = matrix.newMatrixAnalytic(phi_rotate, betta_rotate, phi_mu, betta_mu)
        e_obs_mu = np.dot(A_matrix_analytic, e_obs)  # переход в магнитную СК

        phi_obs, theta_obs = get_angles_from_vector(e_obs_mu)
        # print("theta_obs%d = %f" % (i1, theta_obs))
        for j in range(N_theta_accretion):
            delta_phi_lim = get_lim_for_analytic_integral_phi(theta_range[j], e_obs_mu)
            for i in range(N_phi_accretion):
                if (phi_range[i] >= delta_phi_lim) and (phi_range[i] <= 2 * np.pi - delta_phi_lim):
                    cos_psi_range[j][i] = 1
                else:
                    cos_psi_range[j][i] = -1
        crf[i1] = axes[row_figure, column_figure].contourf(phi_range, theta_range / grad_to_rad,
                                                           cos_psi_range, vmin=-1, vmax=1, cmap=cmap,
                                                           norm=normalizer)
        axes[row_figure, column_figure].set_title(
            "phase = %.2f" % (omega_ns * (t_max_for_cos / (number_of_plots - 1)) * i1 / (2 * np.pi)))

        # axes[row_figure, column_figure].set_theta_zero_location('E', offset=phi_obs/grad_to_rad)
        axes[row_figure, column_figure].set_theta_offset(phi_obs)

        axes[row_figure, column_figure].set_rorigin(-theta_accretion_begin)

        column_figure += 1
        if column_figure == column_number:
            column_figure = 0
            row_figure += 1

    plt.subplots_adjust(hspace=0.5, wspace=0.5)
    cbar = fig.colorbar(im, ax=axes[:, :], shrink=0.7, location='right')
    plt.show()


for i in range(t_max):
    print("%d - " % i, end='')
    print(sum_intense[i], sum_simps_integrate[i])

print("BS total luminosity: ", L_x)
print("Calculated total luminosity: ", calculate_total_luminosity(phi_range, theta_range))
print("difference: Calc/BS = %.5f" % (calculate_total_luminosity(phi_range, theta_range) / L_x))

phi_for_plot = list(omega_ns * i / (2 * np.pi) for i in range(t_max))
fig = plt.figure(figsize=(8, 8))
ax3 = fig.add_subplot(111)
# чтобы максимум был сначала - [position_of_max:], [0:position_of_max]
# ax3.plot(phi_for_plot, np.append(sum_intense[position_of_max:], sum_intense[0:position_of_max]), 'b', label='rectangle')
ax3.plot(phi_for_plot, np.append(sum_simps_integrate[position_of_max:], sum_simps_integrate[0:position_of_max]), 'r',
         label='simps')
ax3.plot(phi_for_plot, np.append(analytic_integral_phi[position_of_max:], analytic_integral_phi[0:position_of_max]),
         'b', marker='*', alpha=0.4,
         label=r"$\phi integrate$")
ax3.set_xlabel('phase')
ax3.set_ylabel("isotropic luminosity, erg/s")
ax3.legend()
# ax3.yscale('log')
plt.yscale('log')
plt.show()

print(M_accretion_rate)
print(H)

file_name = "save%d.txt" % file_count
np.savetxt("phi_for_plot.txt", phi_for_plot)
np.savetxt(file_name, np.append(analytic_integral_phi[position_of_max:], analytic_integral_phi[0:position_of_max]))

row_number = 2
column_number = 3
# plot_map_cos(n_pos, position_of_max, t_max, N_phi_accretion, N_theta_accretion, row_number, column_number)
plot_map_cos_in_range(position_of_max, t_max_for_cos, N_phi_accretion, N_theta_accretion, row_number, column_number)
# plot_map_delta_phi(position_of_max, t_max_for_cos, N_phi_accretion, N_theta_accretion, row_number, column_number)

plot_3D_flag = False
# рисуем 3D

if (plot_3D_flag):
    fig = plt.figure(figsize=(8, 8))
    ax = plt.axes(projection='3d')

    # рисуем звезду
    theta_range = np.arange(0, np.pi, np.pi / N_theta_accretion)
    phi_range = np.arange(0, 2 * np.pi, 2 * np.pi / N_phi_accretion)

    u, v = np.meshgrid(phi_range, theta_range)
    r1 = np.sin(theta_accretion_begin) ** 2
    x = r1 * np.sin(v) * np.cos(u)
    y = r1 * np.sin(v) * np.sin(u)
    z = r1 * np.cos(v)

    ax.plot_surface(x, y, z, color='b', alpha=1)

    # рисуем силовые линии
    theta_range = np.arange(0, np.pi, np.pi / N_theta_accretion)
    phi_range = np.arange(0, 1 / 2 * np.pi, 1 / 2 * np.pi / N_phi_accretion)

    r, p = np.meshgrid(np.sin(theta_range) ** 2, phi_range)
    r1 = r * np.sin(theta_range)
    x = r1 * np.cos(p)
    y = r1 * np.sin(p)
    z = r * np.cos(theta_range)

    ax.plot_wireframe(x, y, z, rstride=4, cstride=4, color="r", alpha=0.2)
    # ax.plot_surface(x, y, z, cmap=plt.cm.YlGnBu_r)

    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])
    ax.set_zlim([-1, 1])
    plt.show()
