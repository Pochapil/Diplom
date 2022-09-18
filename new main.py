import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate

import geometricTask.matrix as matrix
import vectors
import BS_distribution_T_eff as get_T_eff
import config

# const
# альфвеновский радиус и радиус магнитного диполя
R_alfven = (config.mu ** 2 / (2 * config.M_accretion_rate * (2 * config.G * config.M_ns) ** (1 / 2))) ** (2 / 7)
R_e = config.ksi_param * R_alfven  # между 1 и 2 формулой в статье
R_e_outer_surface, R_e_inner_surface = R_e, R_e
# вектор на наблюдателя в системе координат двойной системы
e_obs = np.array([0, np.sin(config.i_angle), np.cos(config.i_angle)])


# формула 2 в статье
def get_delta_distance(theta):
    # R=R_e * sin_theta ** 2
    return R_e * np.sin(theta) ** 3 / (1 + 3 * np.cos(theta) ** 2) ** (1 / 2) * config.dRe_div_Re


# формула 3 в статье
def get_A_normal(theta):
    # a - в азимутальном направлении поток занимает фиксированную долю a полного круга 2πR sinθ
    return 2 * get_delta_distance(theta) * 2 * np.pi * config.a_portion * R_e * np.sin(theta) ** 3


def get_theta_accretion_begin(R_e):
    return np.arcsin((config.R_ns / R_e) ** (1 / 2))


def check_if_intersect(origin_phi, origin_theta, direction_vector, ksi_shock, lim_theta_top, lim_theta_bot, flag):
    # сначала поиск со сферой после - с конусом (а нужно ли со сферой искать ?)
    # все нормирую на радиус НЗ
    # r = origin + t * direction - уравнение луча
    r = R_e * np.sin(origin_theta) ** 2 / config.R_ns

    x_origin = np.sin(origin_theta) * np.cos(origin_phi) * r
    y_origin = np.sin(origin_theta) * np.sin(origin_phi) * r
    z_origin = np.cos(origin_theta) * r
    # origin_point = [x, y, z]

    x_direction = direction_vector[0, 0]
    y_direction = direction_vector[0, 1]
    z_direction = direction_vector[0, 2]

    def find_intersect_solution(a, b, c):
        if b ** 2 - 4 * a * c >= 0:
            t_1 = (-b + (b ** 2 - 4 * a * c) ** (1 / 2)) / (2 * a)
            t_2 = (-b - (b ** 2 - 4 * a * c) ** (1 / 2)) / (2 * a)
            return t_1, t_2
        else:
            return -1, -1

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
    # ограничиваю 2 конусами в зависимости от поверхности (чтобы не закрыть большую часть аппроксимацией)
    #
    lim_theta = lim_theta_top
    if flag:
        lim_theta = lim_theta_bot

    a_cone = (x_direction ** 2 + y_direction ** 2) / (np.tan(lim_theta) ** 2) - z_direction ** 2
    b_cone = 2 * ((x_origin * x_direction + y_origin * y_direction) / (np.tan(lim_theta) ** 2) - z_origin * z_direction)
    c_cone = (x_origin ** 2 + y_origin ** 2) / (np.tan(lim_theta) ** 2) - z_origin ** 2

    # solution for cone
    t_cone = find_intersect_solution(a_cone, b_cone, c_cone)
    # print("t_cone1 = %f,t_cone2 = %f" % (t_cone[0], t_cone[1]))

    for t in t_cone:
        if t > 0:
            intersect_point = np.array([x_origin, y_origin, z_origin]) + t * np.array(
                [x_direction, y_direction, z_direction])
            phi_intersect, theta_intersect = vectors.get_angles_from_vector_one_dimension(intersect_point)
            # для верхнего конуса:
            if (0 < intersect_point[2] < ksi_shock * np.cos(
                    lim_theta_top) and phi_intersect < config.lim_phi_accretion):
                return True
            # для нижнего конуса:
            if -ksi_shock * np.cos(lim_theta_top) < intersect_point[2] < 0:
                if (0 < phi_intersect < config.lim_phi_accretion - np.pi) or (
                        np.pi < phi_intersect < (config.lim_phi_accretion + np.pi)):
                    return True
    return False


class AccretionColumn:

    def __init__(self, R_e_outer_surface, theta_accretion_begin_outer_surface, R_e_inner_surface,
                 theta_accretion_begin_inner_surface, flag):
        self.flag = flag
        self.outer_surface = self.Surface(theta_accretion_begin_outer_surface, R_e_outer_surface, True, self.flag)
        self.inner_surface = self.Surface(theta_accretion_begin_inner_surface, R_e_inner_surface, False, self.flag)

    class Surface:
        def __init__(self, theta_accretion_begin, R_e, flag, column_flag):
            self.R_e = R_e
            self.flag = flag  # True - внешняя поверхность, False - внутренняя

            self.T_eff, self.ksi_shock, self.L_x = get_T_eff.get_Teff_distribution(config.N_theta_accretion, R_e,
                                                                                   get_delta_distance(
                                                                                       theta_accretion_begin),
                                                                                   get_A_normal(theta_accretion_begin))
            phi_delta = 0
            if column_flag:
                theta_accretion_end = np.arcsin((config.R_ns * self.ksi_shock / R_e) ** (1 / 2))
            else:
                theta_accretion_end = np.pi - np.arcsin((config.R_ns * self.ksi_shock / R_e) ** (1 / 2))
                phi_delta = np.pi
            step_phi_accretion = config.lim_phi_accretion / config.N_phi_accretion
            step_theta_accretion = (theta_accretion_end - theta_accretion_begin) / config.N_theta_accretion

            self.theta_range = np.array(
                [theta_accretion_begin + step_theta_accretion * j for j in range(config.N_theta_accretion)])
            self.phi_range = np.array(
                [config.phi_accretion_begin + phi_delta + step_phi_accretion * i for i in
                 range(config.N_phi_accretion)])

            self.cos_psi_range = []  # тут создать матрицу косинусов 1 раз и использовать потом
            self.array_normal = self.create_array_normal(self.phi_range, self.theta_range, self.flag)

        def create_array_normal(self, phi_range, theta_range, flag=True):
            array_normal = []  # матрица нормалей
            coefficient = -1
            if flag:  # True - внешняя поверхность, False - внутренняя
                coefficient = 1
            for i in range(config.N_phi_accretion):
                for j in range(config.N_theta_accretion):
                    array_normal.append(coefficient * matrix.newE_n(phi_range[i], theta_range[j]))
            return array_normal

        def fill_cos_psi_range(self, theta_accretion_begin, theta_accretion_end):
            # sum_intense изотропная светимость ( * 4 pi еще надо)
            # для интеграла по simpson
            cos_psi_range_final = []
            simps_cos = [0] * config.N_theta_accretion  # cos для интеграла по симпсону
            for t in range(config.t_max):
                cos_psi_range = np.empty([config.N_phi_accretion, config.N_theta_accretion])
                # поворот
                phi_mu = config.phi_mu_0 + config.omega_ns * config.grad_to_rad * t
                # расчет матрицы поворота в магнитную СК и вектора на наблюдателя
                A_matrix_analytic = matrix.newMatrixAnalytic(config.phi_rotate, config.betta_rotate, phi_mu,
                                                             config.betta_mu)
                e_obs_mu = np.dot(A_matrix_analytic, e_obs)  # переход в магнитную СК

                # sum_intense изотропная светимость ( * 4 pi еще надо)
                for i in range(config.N_phi_accretion):
                    for j in range(config.N_theta_accretion):
                        # умножать на N_theta
                        cos_psi_range[i, j] = np.dot(e_obs_mu, self.array_normal[i * config.N_theta_accretion + j])
                        if cos_psi_range[i, j] > 0:
                            if check_if_intersect(self.phi_range[i], self.theta_range[j], e_obs_mu, self.ksi_shock,
                                                  theta_accretion_begin, theta_accretion_end, self.flag):
                                cos_psi_range[i, j] = 0
                        else:
                            cos_psi_range[i, j] = 0
                cos_psi_range_final.append(cos_psi_range)
            self.cos_psi_range = cos_psi_range_final

        def calculate_integral_distribution(self):
            # для интеграла по simpson
            sum_simps_integrate = [0] * config.t_max
            simps_integrate_step = [0] * config.N_phi_accretion
            dS_simps = []
            for j in range(config.N_theta_accretion):
                # R=R_e * sin_theta ** 2; R_phi = R * sin_theta
                dl_simps = self.R_e * (3 * np.cos(self.theta_range[j]) ** 2 + 1) ** (1 / 2) * np.sin(
                    self.theta_range[j])
                dphi_simps = self.R_e * np.sin(self.theta_range[j]) ** 3
                dS_simps.append(dphi_simps * dl_simps)  # единичная площадка при интегрировании

            for t in range(config.t_max):
                for i in range(config.N_phi_accretion):
                    simps_integrate_step[i] = np.abs(config.sigmStfBolc * scipy.integrate.simps(
                        self.T_eff ** 4 * np.array(self.cos_psi_range[t][i][:]) * np.array(dS_simps), self.theta_range))
                sum_simps_integrate[t] = scipy.integrate.simps(simps_integrate_step, self.phi_range)

            return sum_simps_integrate


# от поверхности NS - угол при котором радиус = радиусу НЗ
# ----------------- начало инициализации верхней колонки ------------------------
theta_accretion_begin_outer_surface = get_theta_accretion_begin(R_e_outer_surface)
theta_accretion_begin_inner_surface = get_theta_accretion_begin(R_e_inner_surface)

top_column = AccretionColumn(R_e_outer_surface, theta_accretion_begin_outer_surface, R_e_inner_surface,
                             theta_accretion_begin_inner_surface, True)
# ----------------- конец инициализации верхней колонки ------------------------

# ----------------- начало инициализации нижней колонки ------------------------
theta_accretion_begin_outer_surface = np.pi - theta_accretion_begin_outer_surface
theta_accretion_begin_inner_surface = np.pi - theta_accretion_begin_inner_surface

bot_column = AccretionColumn(R_e_outer_surface, theta_accretion_begin_outer_surface, R_e_inner_surface,
                             theta_accretion_begin_inner_surface, False)
# ----------------- конец инициализации нижней колонки ------------------------

# ----------------- углы для нахождения пересечений -------------------------
theta_accretion_begin = top_column.outer_surface.theta_range[0]
theta_accretion_end = top_column.outer_surface.theta_range[-1]
# ---------------------------------------------------------------------------

# ------------------ начало заполнения матриц косинусов ---------------------------
top_column.outer_surface.fill_cos_psi_range(theta_accretion_begin, theta_accretion_end)
top_column.inner_surface.fill_cos_psi_range(theta_accretion_begin, theta_accretion_end)

bot_column.outer_surface.fill_cos_psi_range(theta_accretion_begin, theta_accretion_end)
bot_column.inner_surface.fill_cos_psi_range(theta_accretion_begin, theta_accretion_end)
# ------------------ конец заполнения матриц косинусов ---------------------------

arr_sum_simps_integrate = [0] * 4
i = 0
arr_sum_simps_integrate[i] = top_column.outer_surface.calculate_integral_distribution()
i += 1
arr_sum_simps_integrate[i] = top_column.inner_surface.calculate_integral_distribution()
i += 1
arr_sum_simps_integrate[i] = bot_column.outer_surface.calculate_integral_distribution()
i += 1
arr_sum_simps_integrate[i] = bot_column.inner_surface.calculate_integral_distribution()
i += 1

sum_simps_integrate = np.array(arr_sum_simps_integrate[0])
for i in range(1, 4):
    sum_simps_integrate += np.array(arr_sum_simps_integrate[i])

print('phi_theta_range saved')
file_name = "save_phi_range.txt"
np.savetxt(file_name, top_column.outer_surface.phi_range)
file_name = "save_theta_range.txt"
np.savetxt(file_name, top_column.outer_surface.theta_range)

fig = plt.figure(figsize=(8, 8))
phi_for_plot = list(config.omega_ns * config.grad_to_rad * i / (2 * np.pi) for i in range(config.t_max))
ax3 = fig.add_subplot(111)
ax3.plot(phi_for_plot, arr_sum_simps_integrate[0],
         label='top outer')
ax3.plot(phi_for_plot, arr_sum_simps_integrate[1],
         label='top inner')
ax3.plot(phi_for_plot, arr_sum_simps_integrate[2],
         label='bot outer', marker='*')
ax3.plot(phi_for_plot, arr_sum_simps_integrate[3],
         label='bot inner')
ax3.plot(phi_for_plot, sum_simps_integrate,
         label='sum')
ax3.legend()
# plt.yscale('log')
plt.show()
