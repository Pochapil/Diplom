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
print('R_e = %f' % (R_e / config.R_ns))
R_e_outer_surface, R_e_inner_surface = R_e, R_e
# вектор на наблюдателя в системе координат двойной системы
e_obs = np.array([0, np.sin(config.i_angle), np.cos(config.i_angle)])
file_name_variables = "betta_omega=%d betta_mu=%d a_portion=%f M_rate_c2_Led=%d" \
                      % (config.betta_rotate, config.betta_mu, config.a_portion, config.M_rate_c2_Led)

approx_type = True
approx_method = "cone"
if approx_type:
    approx_method = "dipole"

file_folder = 'figs/'


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


def check_if_intersect(origin_phi, origin_theta, direction_vector, ksi_shock, theta_accretion_begin,
                       theta_accretion_end, top_column_phi_range, bot_column_phi_range, surface_type):
    # сначала поиск со сферой после - с конусом (а нужно ли со сферой искать ?)
    # все нормирую на радиус НЗ
    # r = origin + t * direction - уравнение луча

    r = R_e / config.R_ns * np.sin(origin_theta) ** 2
    # декартовая СК из сферических
    origin_x = np.sin(origin_theta) * np.cos(origin_phi) * r
    origin_y = np.sin(origin_theta) * np.sin(origin_phi) * r
    origin_z = np.cos(origin_theta) * r
    # origin_point = [x, y, z]

    direction_x = direction_vector[0, 0]
    direction_y = direction_vector[0, 1]
    direction_z = direction_vector[0, 2]

    direction_phi, direction_theta = vectors.get_angles_from_vector(direction_vector)

    def find_intersect_solution(a, b, c):
        if b ** 2 - 4 * a * c >= 0:
            t_1 = (-b + (b ** 2 - 4 * a * c) ** (1 / 2)) / (2 * a)
            t_2 = (-b - (b ** 2 - 4 * a * c) ** (1 / 2)) / (2 * a)
            return t_1, t_2
        else:
            return -1, -1

    def intersection_with_sphere():
        # sphere x**2 + y**2 + z**2 == 1
        a_sphere = direction_x ** 2 + direction_y ** 2 + direction_z ** 2
        b_sphere = 2 * (origin_x * direction_x + origin_y * direction_y + origin_z * direction_z)
        c_sphere = origin_x ** 2 + origin_y ** 2 + origin_z ** 2 - 1

        # solution for sphere
        t_sphere = find_intersect_solution(a_sphere, b_sphere, c_sphere)

        # print("t_sphere1 = %f,t_sphere2 = %f" % (t_sphere[0], t_sphere[1]))
        if t_sphere[0] > 0 or t_sphere[1] > 0:
            return True

        return False

    def intersection_with_cone():
        # cone x**2 + y**2 == z**2
        # ограничиваю 2 конусами в зависимости от поверхности (чтобы не закрыть большую часть аппроксимацией)
        #
        lim_theta = theta_accretion_begin
        if surface_type:
            lim_theta = theta_accretion_end

        a_cone = (direction_x ** 2 + direction_y ** 2) / (np.tan(lim_theta) ** 2) - direction_z ** 2
        b_cone = 2 * ((origin_x * direction_x + origin_y * direction_y) / (
                np.tan(lim_theta) ** 2) - origin_z * direction_z)
        c_cone = (origin_x ** 2 + origin_y ** 2) / (np.tan(lim_theta) ** 2) - origin_z ** 2

        # solution for cone
        t_cone = find_intersect_solution(a_cone, b_cone, c_cone)
        # print("t_cone1 = %f,t_cone2 = %f" % (t_cone[0], t_cone[1]))

        for t in t_cone:
            if t > 0:
                intersect_point = np.array([origin_x, origin_y, origin_z]) + t * np.array(
                    [direction_x, direction_y, direction_z])
                intersect_phi, intersect_theta = vectors.get_angles_from_vector_one_dimension(intersect_point)

                # для верхнего конуса:
                intersect_z_correct = 0 < intersect_point[2] < ksi_shock * np.cos(theta_accretion_end)
                intersect_phi_correct = top_column_phi_range[0] < intersect_phi < top_column_phi_range[-1]

                if intersect_z_correct and intersect_phi_correct:
                    return True

                # для нижнего конуса:
                intersect_z_correct = -ksi_shock * np.cos(theta_accretion_end) < intersect_point[2] < 0
                # условие - так как углы из метода get_angles_from_vector_one_dimension от 0 до 2 pi поэтому
                # нужно учесть углы которые превышают 2 pi в 1 скобке условия
                intersect_phi_correct = (0 < intersect_phi < bot_column_phi_range[-1] - 2 * np.pi) or (
                        bot_column_phi_range[0] < intersect_phi < bot_column_phi_range[-1])

                if intersect_z_correct and intersect_phi_correct:
                    return True
        return False

    def intersection_with_dipole_lines():
        '''
        есть аналитическое уравнение для полинома дипольной линии 5 степени
        находим уравнение в сферических координатах.

        достаем корни
        ищем положительные
        находим пересечение с колонками
            если пересекается то истина
            если нет то ложь
        '''
        # для 0 угла наблюдателя по фи

        phi_delta = origin_phi - direction_phi
        eta = np.sin(direction_theta) / np.sin(origin_theta)
        cos_alpha = np.sin(origin_theta) * np.cos(phi_delta) * np.sin(direction_theta) + np.cos(origin_theta) * np.cos(
            direction_theta)

        c_x_5 = 1
        c_x_4 = 6 * cos_alpha
        c_x_3 = 3 + 12 * cos_alpha ** 2 - eta ** 4
        c_x_2 = 12 * cos_alpha + 8 * cos_alpha ** 3 - 4 * np.cos(phi_delta) * eta ** 3
        c_x_1 = 3 + 12 * cos_alpha ** 2 - 2 * eta ** 2 - 4 * np.cos(phi_delta) ** 2 * eta ** 2
        c_x_0 = 6 * cos_alpha - 4 * np.cos(phi_delta) * eta

        coefficients = [c_x_5, c_x_4, c_x_3, c_x_2, c_x_1, c_x_0]
        solutions = np.roots(coefficients)

        for solution in solutions:
            if solution.real > 0 and solution.imag == 0:
                direction_t = solution.real * r
                intersect_point = np.array([origin_x, origin_y, origin_z]) + direction_t * np.array(
                    [direction_x, direction_y, direction_z])
                intersect_phi, intersect_theta = vectors.get_angles_from_vector_one_dimension(intersect_point)
                intersect_r = (intersect_point[0] ** 2 + intersect_point[1] ** 2 + intersect_point[2] ** 2) ** (1 / 2)

                intersect_r_correct = intersect_r < ksi_shock

                # для верхнего конуса:
                intersect_z_correct = 0 < intersect_point[2] < ksi_shock * np.cos(theta_accretion_end)
                intersect_phi_correct = top_column_phi_range[0] < intersect_phi < top_column_phi_range[-1]
                if intersect_z_correct and intersect_phi_correct and intersect_r_correct:
                    return True

                # для нижнего конуса:
                intersect_z_correct = - (ksi_shock * np.cos(theta_accretion_end)) < intersect_point[2] < 0
                # условие - так как углы из метода get_angles_from_vector_one_dimension от 0 до 2 pi поэтому
                # нужно учесть углы которые превышают 2 pi в 1 скобке условия
                intersect_phi_correct = (0 < intersect_phi < bot_column_phi_range[-1] - 2 * np.pi) or (
                        bot_column_phi_range[0] < intersect_phi < bot_column_phi_range[-1])
                if intersect_z_correct and intersect_phi_correct and intersect_r_correct:
                    return True

        return False

    if approx_type:
        return (intersection_with_sphere() or intersection_with_dipole_lines())
    else:
        return (intersection_with_sphere() or intersection_with_cone())


def get_pulsed_fraction(arr):
    min_value = min(arr)
    max_value = max(arr)
    return (max_value - min_value) / (max_value + min_value)


class AccretionColumn:

    def __init__(self, R_e_outer_surface, theta_accretion_begin_outer_surface, R_e_inner_surface,
                 theta_accretion_begin_inner_surface, column_type):
        self.column_type = column_type
        self.outer_surface = self.Surface(theta_accretion_begin_outer_surface, R_e_outer_surface, True,
                                          self.column_type)
        self.inner_surface = self.Surface(theta_accretion_begin_inner_surface, R_e_inner_surface, False,
                                          self.column_type)

    class Surface:
        def __init__(self, theta_accretion_begin, R_e, surface_type, column_type):
            self.R_e = R_e
            self.surface_type = surface_type  # True - внешняя поверхность, False - внутренняя

            self.T_eff, self.ksi_shock, self.L_x = get_T_eff.get_Teff_distribution(config.N_theta_accretion, R_e,
                                                                                   get_delta_distance(
                                                                                       theta_accretion_begin),
                                                                                   get_A_normal(theta_accretion_begin))
            phi_delta = 0
            theta_accretion_end = np.arcsin((config.R_ns * self.ksi_shock / R_e) ** (1 / 2))
            if not column_type:
                theta_accretion_end = np.pi - theta_accretion_end
                phi_delta = np.pi
            step_phi_accretion = config.lim_phi_accretion / config.N_phi_accretion
            step_theta_accretion = (theta_accretion_end - theta_accretion_begin) / config.N_theta_accretion

            self.theta_range = np.array(
                [theta_accretion_begin + step_theta_accretion * j for j in range(config.N_theta_accretion)])
            self.phi_range = np.array(
                [config.phi_accretion_begin + phi_delta + step_phi_accretion * i for i in
                 range(config.N_phi_accretion)])

            self.cos_psi_range = []  # тут создать матрицу косинусов 1 раз и использовать потом
            self.array_normal = self.create_array_normal(self.phi_range, self.theta_range, self.surface_type)

            self.correct_T_eff()

        def correct_T_eff(self):
            # так как нахожу распределение Teff по отрезку кси, то нужно привести к виду по сетке theta !!
            # нужно сместить для увеличения точности
            ksi_stop = 1.
            ksi_inc = - (self.ksi_shock - ksi_stop) / config.N_theta_accretion
            ksi_range = np.arange(self.ksi_shock, ksi_stop, ksi_inc)
            ksi_bs = ksi_range[::-1]
            i = 0  # left_border
            true_T_eff = []
            for theta in self.theta_range:
                while ksi_bs[i + 1] < self.R_e / config.R_ns * np.sin(theta) ** 2 and (
                        i < config.N_theta_accretion - 2):
                    i += 1
                true_T_eff.append(self.T_eff[i])
            self.T_eff = np.array(true_T_eff)

        def create_array_normal(self, phi_range, theta_range, surface_type=True):
            array_normal = []  # матрица нормалей
            coefficient = -1
            if surface_type:  # True - внешняя поверхность, False - внутренняя
                coefficient = 1
            for i in range(config.N_phi_accretion):
                for j in range(config.N_theta_accretion):
                    array_normal.append(coefficient * matrix.newE_n(phi_range[i], theta_range[j]))
            return array_normal

        def fill_cos_psi_range(self, theta_accretion_begin, theta_accretion_end, top_column_phi_range,
                               bot_column_phi_range):
            # sum_intense изотропная светимость ( * 4 pi еще надо)
            # для интеграла по simpson
            cos_psi_range_final = []
            for t in range(config.t_max):
                cos_psi_range = np.empty([config.N_phi_accretion, config.N_theta_accretion])
                # поворот
                phi_mu = config.phi_mu_0 + config.omega_ns * config.grad_to_rad * t
                # расчет матрицы поворота в магнитную СК и вектора на наблюдателя
                A_matrix_analytic = matrix.newMatrixAnalytic(config.phi_rotate, config.betta_rotate, phi_mu,
                                                             config.betta_mu)
                e_obs_mu = np.dot(A_matrix_analytic, e_obs)  # переход в магнитную СК
                for i in range(config.N_phi_accretion):
                    for j in range(config.N_theta_accretion):
                        # умножать на N_theta
                        cos_psi_range[i, j] = np.dot(e_obs_mu, self.array_normal[i * config.N_theta_accretion + j])
                        if cos_psi_range[i, j] > 0:
                            if check_if_intersect(self.phi_range[i], self.theta_range[j], e_obs_mu, self.ksi_shock,
                                                  theta_accretion_begin, theta_accretion_end, top_column_phi_range,
                                                  bot_column_phi_range, self.surface_type):
                                cos_psi_range[i, j] = 0
                        else:
                            cos_psi_range[i, j] = 0
                cos_psi_range_final.append(cos_psi_range)
            self.cos_psi_range = cos_psi_range_final

        def calculate_integral_distribution(self):
            # для интеграла по simpson
            dS_simps = []  # единичная площадка при интегрировании
            for j in range(config.N_theta_accretion):
                # R=R_e * sin_theta ** 2; R_phi = R * sin_theta
                dl_simps = self.R_e * (3 * np.cos(self.theta_range[j]) ** 2 + 1) ** (1 / 2) * np.sin(
                    self.theta_range[j])
                dphi_simps = self.R_e * np.sin(self.theta_range[j]) ** 3
                dS_simps.append(dphi_simps * dl_simps)  # единичная площадка при интегрировании

            sum_simps_integrate = [0] * config.t_max
            simps_integrate_step = [0] * config.N_phi_accretion
            for t in range(config.t_max):
                for i in range(config.N_phi_accretion):
                    # /pi т.к считаем мощность, которая приходит на наблюдателя, а не всю светимость (пи входит в сигма)
                    simps_integrate_step[i] = np.abs(config.sigmStfBolc / np.pi * scipy.integrate.simps(
                        self.T_eff ** 4 * np.array(self.cos_psi_range[t][i][:]) * np.array(dS_simps), self.theta_range))
                sum_simps_integrate[t] = scipy.integrate.simps(simps_integrate_step, self.phi_range)

            return sum_simps_integrate

        def calculate_integral_distribution_in_range(self, energy_bot, energy_top):
            # КэВ
            '''
            на каждом слое по тета берем в данном диапазоне интеграл по формуле планка
            переходя на след слой меняем температуру
            суммируем
            конец
            
            внешний цикл - по длине волны
            внутри - идем в цикле по углам тета
            '''

            # нужно ли на пи?
            def plank_energy_on_wavelength(wavelength, T):
                return 2 * config.h_plank_ergs * config.c ** 2 / wavelength ** 5 \
                       * 1 / (np.e ** (config.h_plank_ergs * config.c / (wavelength * config.k_bolc * T)) - 1)

            def plank_energy_on_frequency(frequency, T):
                return 2 * config.h_plank_ergs * frequency ** 3 / config.c ** 2 \
                       * 1 / (np.e ** (config.h_plank_ergs * frequency / (config.k_bolc * T)) - 1)

            coefficient = 1000  # КэВ а не эВ
            frequency_top = coefficient * energy_top / config.h_plank_evs  # E = h f
            frequency_bot = coefficient * energy_bot / config.h_plank_evs  # E = h f
            frequency_step = (frequency_top - frequency_bot) / config.N_frequency_range
            frequency_range = [frequency_bot + frequency_step * i for i in range(config.N_frequency_range)]

            dS_simps = []  # единичная площадка при интегрировании
            for j in range(config.N_theta_accretion):
                # R=R_e * sin_theta ** 2; R_phi = R * sin_theta
                dl_simps = self.R_e * (3 * np.cos(self.theta_range[j]) ** 2 + 1) ** (1 / 2) * np.sin(
                    self.theta_range[j])
                dphi_simps = self.R_e * np.sin(self.theta_range[j]) ** 3
                dS_simps.append(dphi_simps * dl_simps)  # единичная площадка при интегрировании

            plank_func = [0] * config.N_theta_accretion
            plank_func_step = [0] * config.N_frequency_range
            for theta_index in range(config.N_theta_accretion):
                for frequency_index in range(config.N_frequency_range):
                    plank_func_step[frequency_index] = plank_energy_on_frequency(frequency_range[frequency_index],
                                                                                 self.T_eff[theta_index])
                plank_func[theta_index] = scipy.integrate.simps(plank_func_step, frequency_range)

            integrate_step = [0] * config.N_phi_accretion
            integrate_sum = [0] * config.t_max
            for rotation_index in range(config.t_max):
                for phi_index in range(config.N_phi_accretion):
                    integrate_step[phi_index] = np.abs(scipy.integrate.simps(
                        plank_func * np.array(dS_simps) * np.array(
                            self.cos_psi_range[rotation_index][phi_index][:]), self.theta_range))
                integrate_sum[rotation_index] = np.abs(scipy.integrate.simps(integrate_step, self.phi_range))

            return integrate_sum


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

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
ax.plot(top_column.outer_surface.theta_range, top_column.outer_surface.T_eff)
ax.legend()
plt.show()


bot_column = AccretionColumn(R_e_outer_surface, theta_accretion_begin_outer_surface, R_e_inner_surface,
                             theta_accretion_begin_inner_surface, False)
# ----------------- конец инициализации нижней колонки ------------------------

print('phi_theta_range saved')
file_name = "save_phi_range.txt"
full_file_name = file_folder + file_name
np.savetxt(full_file_name, top_column.outer_surface.phi_range)
file_name = "save_theta_range.txt"
full_file_name = file_folder + file_name
np.savetxt(full_file_name, top_column.outer_surface.theta_range)

# print('T_eff:')
# print(top_column.outer_surface.T_eff)

# ----------------- углы для нахождения пересечений -------------------------
theta_accretion_begin = top_column.outer_surface.theta_range[0]
theta_accretion_end = top_column.outer_surface.theta_range[-1]
# ---------------------------------------------------------------------------

# ------------------ начало заполнения матриц косинусов ---------------------------
top_column.outer_surface.fill_cos_psi_range(theta_accretion_begin, theta_accretion_end,
                                            top_column.outer_surface.phi_range, bot_column.outer_surface.phi_range)
top_column.inner_surface.fill_cos_psi_range(theta_accretion_begin, theta_accretion_end,
                                            top_column.outer_surface.phi_range, bot_column.outer_surface.phi_range)

bot_column.outer_surface.fill_cos_psi_range(theta_accretion_begin, theta_accretion_end,
                                            top_column.outer_surface.phi_range, bot_column.outer_surface.phi_range)
bot_column.inner_surface.fill_cos_psi_range(theta_accretion_begin, theta_accretion_end,
                                            top_column.outer_surface.phi_range, bot_column.outer_surface.phi_range)
# ------------------ конец заполнения матриц косинусов ---------------------------

arr_simps_integrate = [0] * 4
i = 0
arr_simps_integrate[i] = top_column.outer_surface.calculate_integral_distribution()
# file_name = "%s %s %d.txt" % (file_name_variables, approx_method, i)
# np.savetxt(file_name, arr_simps_integrate[i])
i += 1
arr_simps_integrate[i] = top_column.inner_surface.calculate_integral_distribution()
# file_name = "%s %s %d.txt" % (file_name_variables, approx_method, i)
# np.savetxt(file_name, arr_simps_integrate[i])
i += 1
arr_simps_integrate[i] = bot_column.outer_surface.calculate_integral_distribution()
# file_name = "%s %s %d.txt" % (file_name_variables, approx_method, i)
# np.savetxt(file_name, arr_simps_integrate[i])
i += 1
arr_simps_integrate[i] = bot_column.inner_surface.calculate_integral_distribution()
# file_name = "%s %s %d.txt" % (file_name_variables, approx_method, i)
# np.savetxt(file_name, arr_simps_integrate[i])
i += 1

print('ksi_shock = %f' % bot_column.outer_surface.ksi_shock)

sum_simps_integrate = np.array(arr_simps_integrate[0])
for i in range(1, 4):
    sum_simps_integrate += np.array(arr_simps_integrate[i])

fig = plt.figure(figsize=(8, 8))
phi_for_plot = list(config.omega_ns * config.grad_to_rad * i / (2 * np.pi) for i in range(config.t_max_for_plot))

append_index = config.t_max_for_plot - config.t_max
for i in range(4):
    arr_simps_integrate[i] = np.append(arr_simps_integrate[i], arr_simps_integrate[i][0:append_index])
sum_simps_integrate = np.append(sum_simps_integrate, sum_simps_integrate[0:append_index])

ax = fig.add_subplot(111)
ax.plot(phi_for_plot, arr_simps_integrate[0],
        label='top outer')
ax.plot(phi_for_plot, arr_simps_integrate[1],
        label='top inner')
ax.plot(phi_for_plot, arr_simps_integrate[2],
        label='bot outer', marker='*')
ax.plot(phi_for_plot, arr_simps_integrate[3],
        label='bot inner')
ax.plot(phi_for_plot, sum_simps_integrate,
        label='sum')
ax.legend()
fig.suptitle('total luminosity of surfaces', fontsize=14)
# plt.yscale('log')
file_name = 'total_luminosity_of_surfaces.png'
full_file_name = file_folder + file_name
fig.savefig(full_file_name + 'total_luminosity_of_surfaces.png', dpi=fig.dpi)
plt.show()

observer_theta = [0] * config.t_max_for_plot
observer_phi = [0] * config.t_max_for_plot

for t in range(config.t_max_for_plot):
    phi_mu = config.phi_mu_0 + config.omega_ns * config.grad_to_rad * t
    # расчет матрицы поворота в магнитную СК и вектора на наблюдателя
    A_matrix_analytic = matrix.newMatrixAnalytic(config.phi_rotate, config.betta_rotate, phi_mu,
                                                 config.betta_mu)
    e_obs_mu = np.dot(A_matrix_analytic, e_obs)  # переход в магнитную СК
    observer_phi[t], observer_theta[t] = vectors.get_angles_from_vector(e_obs_mu)

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
ax.plot(phi_for_plot, observer_theta, label=r'$\theta_{observer}$')
ax.plot(phi_for_plot, observer_phi, label=r'$\phi_{observer}$')
ax.legend()
fig.suptitle('Observer angles', fontsize=14)
file_name = 'Observer_angles.png'
full_file_name = file_folder + file_name
fig.savefig(full_file_name + 'Observer_angles.png', dpi=fig.dpi)
plt.show()

while True:
    energy_bot = float(input('введите нижний предел в КэВ: '))
    energy_top = float(input('введите верхний предел в КэВ: '))

    arr_simps_integrate = [0] * 4
    i = 0
    arr_simps_integrate[i] = top_column.outer_surface.calculate_integral_distribution_in_range(energy_bot, energy_top)
    i += 1
    arr_simps_integrate[i] = top_column.inner_surface.calculate_integral_distribution_in_range(energy_bot, energy_top)
    i += 1
    arr_simps_integrate[i] = bot_column.outer_surface.calculate_integral_distribution_in_range(energy_bot, energy_top)
    i += 1
    arr_simps_integrate[i] = bot_column.inner_surface.calculate_integral_distribution_in_range(energy_bot, energy_top)
    i += 1

    sum_simps_integrate = np.array(arr_simps_integrate[0])
    for i in range(1, 4):
        sum_simps_integrate += np.array(arr_simps_integrate[i])

    PF = get_pulsed_fraction(sum_simps_integrate)

    fig = plt.figure(figsize=(8, 8))
    phi_for_plot = list(config.omega_ns * config.grad_to_rad * i / (2 * np.pi) for i in range(config.t_max_for_plot))

    append_index = config.t_max_for_plot - config.t_max
    for i in range(4):
        arr_simps_integrate[i] = np.append(arr_simps_integrate[i], arr_simps_integrate[i][0:append_index])
    sum_simps_integrate = np.append(sum_simps_integrate, sum_simps_integrate[0:append_index])

    ax = fig.add_subplot(111)
    # ax.plot(phi_for_plot, arr_simps_integrate[0],
    #         label='top outer')
    # ax.plot(phi_for_plot, arr_simps_integrate[1],
    #         label='top inner')
    # ax.plot(phi_for_plot, arr_simps_integrate[2],
    #         label='bot outer', marker='*')
    # ax.plot(phi_for_plot, arr_simps_integrate[3],
    #         label='bot inner')
    ax.plot(phi_for_plot, sum_simps_integrate,
            label='sum')
    ax.legend()
    file_name = "sum_of_luminosity_in_range_%0.2f_-_%0.2f_KeV_of_surfaces.txt" % (energy_bot, energy_top)
    full_file_name = file_folder + file_name
    np.savetxt(full_file_name, sum_simps_integrate)

    fig_title = 'luminosity in range %0.2f - %0.2f KeV of surfaces, PF = %0.3f' % (energy_bot, energy_top, PF)
    file_name = 'luminosity_in_range%0.2f_-_%0.2f_KeV_of_surfaces.png' % (energy_bot, energy_top)
    full_file_name = file_folder + file_name
    fig.suptitle(fig_title, fontsize=14)
    fig.savefig(full_file_name, dpi=fig.dpi)
    # plt.yscale('log')
    plt.show()
