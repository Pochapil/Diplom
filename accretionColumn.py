import numpy as np
import scipy.integrate

import accretionColumnService
import config
import BS_distribution_T_eff as get_T_eff
import geometricTask.matrix as matrix

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
                                                                                   accretionColumnService.get_delta_distance(
                                                                                       theta_accretion_begin, self.R_e),
                                                                                   accretionColumnService.get_A_normal(
                                                                                       theta_accretion_begin, self.R_e))
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
                            if accretionColumnService.check_if_intersect(self.phi_range[i], self.theta_range[j], e_obs_mu, self.ksi_shock,
                                                  theta_accretion_begin, theta_accretion_end, top_column_phi_range,
                                                  bot_column_phi_range, self.surface_type, self.R_e):
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
