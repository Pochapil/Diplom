import numpy as np
import scipy.integrate
from scipy import interpolate

import accretionColumnService
import config
import BS_distribution_T_eff as get_T_eff
import geometricTask.matrix as matrix


class AccretionColumn:

    def __init__(self, R_e_outer_surface, theta_accretion_begin_outer_surface, R_e_inner_surface,
                 theta_accretion_begin_inner_surface, column_type):
        self.column_type = column_type  # True - top, False - bot
        self.outer_surface = self.Surface(theta_accretion_begin_outer_surface, R_e_outer_surface, True,
                                          self.column_type)
        self.inner_surface = self.Surface(theta_accretion_begin_inner_surface, R_e_inner_surface, False,
                                          self.column_type)

    class Surface:
        def __init__(self, theta_accretion_begin, R_e, surface_type, column_type):
            self.R_e = R_e
            self.surface_type = surface_type  # True - внешняя поверхность, False - внутренняя

            self.T_eff, self.ksi_shock, self.L_x = get_T_eff.get_Teff_distribution(config.N_theta_accretion, R_e,
                                                                                   accretionColumnService.get_delta_distance(
                                                                                       theta_accretion_begin, self.R_e),
                                                                                   accretionColumnService.get_A_normal(
                                                                                       theta_accretion_begin, self.R_e))
            phi_delta = 0
            theta_accretion_end = np.arcsin((config.R_ns * self.ksi_shock / R_e) ** (1 / 2))
            if not column_type:  # column_type: True - top, False - bot
                # для нижней сместить углы
                theta_accretion_end = np.pi - theta_accretion_end
                phi_delta = np.pi
            step_phi_accretion = config.lim_phi_accretion / (config.N_phi_accretion - 1)
            step_theta_accretion = (theta_accretion_end - theta_accretion_begin) / (config.N_theta_accretion - 1)

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
            # нужно проинтерполировать внутри отрезка

            # создание сетки по которой было найдено распределение T_eff
            ksi_stop = 1
            ksi_inc = (self.ksi_shock - ksi_stop) / (config.N_theta_accretion - 1)
            ksi_range = np.array([ksi_stop + ksi_inc * i for i in range(config.N_theta_accretion)])

            # интерполяция
            x = ksi_range
            y = self.T_eff
            f = interpolate.interp1d(x, y, kind='cubic')

            # новый диапазон по х - начинаем с 1 индекса так как 0 элемент совпадает - x[0] = 1
            x_new = self.R_e / config.R_ns * np.sin(self.theta_range[1:-1]) ** 2
            y_new = f(x_new)

            # выравнивание
            for i in range(0, (len(y_new))):
                self.T_eff[i + 1] = y_new[i]

        def create_ds_for_integral(self):
            # для интеграла по simpson
            dS_simps = []  # единичная площадка при интегрировании
            for j in range(config.N_theta_accretion):
                # R=R_e * sin_theta ** 2; R_phi = R * sin_theta
                dl_simps = self.R_e * (3 * np.cos(self.theta_range[j]) ** 2 + 1) ** (1 / 2) * np.sin(
                    self.theta_range[j])
                dphi_simps = self.R_e * np.sin(self.theta_range[j]) ** 3
                dS_simps.append(dphi_simps * dl_simps)  # единичная площадка при интегрировании
            return dS_simps

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
                               bot_column_phi_range, e_obs):
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
                            if accretionColumnService.check_if_intersect(self.phi_range[i], self.theta_range[j],
                                                                         e_obs_mu, self.ksi_shock,
                                                                         theta_accretion_begin, theta_accretion_end,
                                                                         top_column_phi_range,
                                                                         bot_column_phi_range, self.surface_type,
                                                                         self.R_e):
                                cos_psi_range[i, j] = 0
                        else:
                            cos_psi_range[i, j] = 0
                cos_psi_range_final.append(cos_psi_range)
            self.cos_psi_range = cos_psi_range_final

        def async_fill_cos_psi_range(self, t_index, theta_accretion_begin, theta_accretion_end, top_column_phi_range,
                                     bot_column_phi_range, e_obs):
            # распараллелил fill_cos_psi_range
            cos_psi_range = np.empty([config.N_phi_accretion, config.N_theta_accretion])
            # поворот
            phi_mu = config.phi_mu_0 + config.omega_ns * config.grad_to_rad * t_index
            # расчет матрицы поворота в магнитную СК и вектора на наблюдателя
            A_matrix_analytic = matrix.newMatrixAnalytic(config.phi_rotate, config.betta_rotate, phi_mu,
                                                         config.betta_mu)
            e_obs_mu = np.dot(A_matrix_analytic, e_obs)  # переход в магнитную СК
            for i in range(config.N_phi_accretion):
                for j in range(config.N_theta_accretion):
                    # умножать на N_theta
                    cos_psi_range[i, j] = np.dot(e_obs_mu, self.array_normal[i * config.N_theta_accretion + j])
                    if cos_psi_range[i, j] > 0:
                        if accretionColumnService.check_if_intersect(self.phi_range[i], self.theta_range[j],
                                                                     e_obs_mu, self.ksi_shock,
                                                                     theta_accretion_begin, theta_accretion_end,
                                                                     top_column_phi_range,
                                                                     bot_column_phi_range, self.surface_type,
                                                                     self.R_e):
                            cos_psi_range[i, j] = 0
                    else:
                        cos_psi_range[i, j] = 0

            return cos_psi_range

        def calculate_integral_distribution(self):
            # luminosity
            # для интеграла по simpson
            dS_simps = self.create_ds_for_integral()

            sum_simps_integrate = [0] * config.t_max
            simps_integrate_step = [0] * config.N_phi_accretion
            for t in range(config.t_max):
                for i in range(config.N_phi_accretion):
                    # /pi т.к считаем мощность, которая приходит на наблюдателя, а не всю светимость (пи входит в сигма)
                    simps_integrate_step[i] = np.abs(config.sigmStfBolc / np.pi * scipy.integrate.simps(
                        self.T_eff ** 4 * np.array(self.cos_psi_range[t][i][:]) * np.array(dS_simps), self.theta_range))
                sum_simps_integrate[t] = scipy.integrate.simps(simps_integrate_step, self.phi_range)

            return sum_simps_integrate

        def calculate_total_luminosity(self):
            dS_simps = self.create_ds_for_integral()

            total_luminosity_step = [0] * config.N_phi_accretion
            for i in range(config.N_phi_accretion):
                total_luminosity_step[i] = config.sigmStfBolc * scipy.integrate.simps(self.T_eff ** 4 * dS_simps,
                                                                                      self.theta_range)
            total_luminosity = scipy.integrate.simps(total_luminosity_step, self.phi_range)
            return total_luminosity

        def async_calculate_integral_distribution(self, t_index):
            # luminosity
            # для интеграла по simpson
            dS_simps = self.create_ds_for_integral()
            simps_integrate_step = [0] * config.N_phi_accretion

            for i in range(config.N_phi_accretion):
                # /pi т.к считаем мощность, которая приходит на наблюдателя, а не всю светимость (пи входит в сигма)
                simps_integrate_step[i] = np.abs(config.sigmStfBolc / np.pi * scipy.integrate.simps(
                    self.T_eff ** 4 * np.array(self.cos_psi_range[t_index][i][:]) * np.array(dS_simps),
                    self.theta_range))
            sum_simps_integrate = scipy.integrate.simps(simps_integrate_step, self.phi_range)

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
                # erg / s / sr / cm**2 / hz
                return 2 * config.h_plank_ergs * frequency ** 3 / config.c ** 2 \
                       * 1 / (np.e ** (config.h_plank_ergs * frequency / (config.k_bolc * T)) - 1)

            coefficient = 1000  # КэВ а не эВ
            frequency_top = coefficient * energy_top / config.h_plank_evs  # E = h f
            frequency_bot = coefficient * energy_bot / config.h_plank_evs  # E = h f
            frequency_step = (frequency_top - frequency_bot) / config.N_frequency_range
            frequency_range = [frequency_bot + frequency_step * i for i in range(config.N_frequency_range)]

            dS_simps = self.create_ds_for_integral()

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

        def calculate_L_nu_on_energy(self, energy):
            # КэВ
            def plank_energy_on_frequency(frequency, T):
                return 2 * config.h_plank_ergs * frequency ** 3 / config.c ** 2 \
                       * 1 / (np.e ** (config.h_plank_ergs * frequency / (config.k_bolc * T)) - 1)

            coefficient = 1000  # КэВ а не эВ
            frequency = coefficient * energy / config.h_plank_evs  # E = h f

            dS_simps = self.create_ds_for_integral()

            plank_func = plank_energy_on_frequency(frequency, self.T_eff)

            integrate_step = [0] * config.N_phi_accretion
            integrate_sum = [0] * config.t_max
            for rotation_index in range(config.t_max):
                for phi_index in range(config.N_phi_accretion):
                    integrate_step[phi_index] = np.abs(scipy.integrate.simps(
                        plank_func * np.array(dS_simps) * np.array(
                            self.cos_psi_range[rotation_index][phi_index][:]), self.theta_range))
                integrate_sum[rotation_index] = np.abs(scipy.integrate.simps(integrate_step, self.phi_range))

            return integrate_sum

        def calculate_nu_L_nu_on_energy(self, energy):
            # КэВ
            coefficient = 1000  # КэВ а не эВ
            frequency = coefficient * energy / config.h_plank_evs  # E = h f
            return np.array(self.calculate_L_nu_on_energy(energy)) * frequency
