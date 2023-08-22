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

            self.T_eff, self.ksi_shock, self.L_x, self.beta = get_T_eff.get_Teff_distribution(R_e,
                                                                                              accretionColumnService.get_delta_distance(
                                                                                                  theta_accretion_begin,
                                                                                                  self.R_e),
                                                                                              accretionColumnService.get_A_normal(
                                                                                                  theta_accretion_begin,
                                                                                                  self.R_e))
            phi_delta = 0

            if (config.R_ns * self.ksi_shock / R_e) > 1:
                # есть набор параметров при которых модель не работает и ударная волна дальше магнитосферы, берем 90
                theta_accretion_end = np.pi / 2
            else:
                # из усл силовой линии МП : r = R_e sin**2; end: ksi_shock = R_e sin**2
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
            if self.ksi_shock > self.R_e / config.R_ns:
                ksi_inc = (self.R_e / config.R_ns - ksi_stop) / (config.N_theta_accretion - 1)
            else:
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
                            # проверка на пересечения
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
                                     bot_column_phi_range, e_obs, updated_betta_mu):
            # распараллелил fill_cos_psi_range
            cos_psi_range = np.empty([config.N_phi_accretion, config.N_theta_accretion])
            # поворот
            phi_mu = config.phi_mu_0 + config.omega_ns * config.grad_to_rad * t_index
            # расчет матрицы поворота в магнитную СК и вектора на наблюдателя
            A_matrix_analytic = matrix.newMatrixAnalytic(config.phi_rotate, config.betta_rotate, phi_mu,
                                                         updated_betta_mu)
            # print('config.betta_mu in fill cos = %f' % updated_betta_mu)
            e_obs_mu = np.dot(A_matrix_analytic, e_obs)  # переход в магнитную СК
            for i in range(config.N_phi_accretion):
                for j in range(config.N_theta_accretion):
                    # умножать на N_theta
                    cos_psi_range[i, j] = np.dot(e_obs_mu, self.array_normal[i * config.N_theta_accretion + j])
                    if cos_psi_range[i, j] > 0:
                        # проверка на пересечения
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
                    # 4 * - для L_iso
                    simps_integrate_step[i] = np.abs(4 * config.sigmStfBolc * scipy.integrate.simps(
                        self.T_eff ** 4 * np.array(self.cos_psi_range[t][i][:]) * np.array(dS_simps), self.theta_range))
                sum_simps_integrate[t] = scipy.integrate.simps(simps_integrate_step, self.phi_range)

            return sum_simps_integrate

        def async_calculate_integral_distribution(self, t_index):
            # luminosity
            # для интеграла по simpson
            dS_simps = self.create_ds_for_integral()
            simps_integrate_step = [0] * config.N_phi_accretion

            for i in range(config.N_phi_accretion):
                # /pi т.к считаем мощность, которая приходит на наблюдателя, а не всю светимость (пи входит в сигма)
                # 4 * - для L_iso
                simps_integrate_step[i] = np.abs(config.sigmStfBolc * scipy.integrate.simps(
                    self.T_eff ** 4 * np.array(self.cos_psi_range[t_index][i][:]) * np.array(dS_simps),
                    self.theta_range))
            sum_simps_integrate = scipy.integrate.simps(simps_integrate_step, self.phi_range)

            return sum_simps_integrate

        def calculate_total_luminosity(self):
            dS_simps = self.create_ds_for_integral()

            total_luminosity_step = [0] * config.N_phi_accretion
            for i in range(config.N_phi_accretion):
                total_luminosity_step[i] = config.sigmStfBolc * scipy.integrate.simps(self.T_eff ** 4 * dS_simps,
                                                                                      self.theta_range)
            total_luminosity = scipy.integrate.simps(total_luminosity_step, self.phi_range)

            return total_luminosity

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

            frequency_top = accretionColumnService.get_frequency_from_energy(energy_top)
            frequency_bot = accretionColumnService.get_frequency_from_energy(energy_bot)
            frequency_range = np.linspace(frequency_bot, frequency_top, config.N_frequency_range)
            dS_simps = self.create_ds_for_integral()

            plank_func = [0] * config.N_theta_accretion
            plank_func_step = [0] * config.N_frequency_range
            for theta_index in range(config.N_theta_accretion):
                for frequency_index in range(config.N_frequency_range):
                    plank_func_step[frequency_index] = accretionColumnService.plank_energy_on_frequency(
                        frequency_range[frequency_index],
                        self.T_eff[theta_index])
                plank_func[theta_index] = scipy.integrate.simps(plank_func_step, frequency_range)

            integrate_step = [0] * config.N_phi_accretion
            integrate_sum = [0] * config.t_max
            for rotation_index in range(config.t_max):
                for phi_index in range(config.N_phi_accretion):
                    # 4 * np.pi для L_iso
                    integrate_step[phi_index] = 4 * np.pi * np.abs(scipy.integrate.simps(
                        plank_func * np.array(dS_simps) * np.array(
                            self.cos_psi_range[rotation_index][phi_index][:]), self.theta_range))
                integrate_sum[rotation_index] = np.abs(scipy.integrate.simps(integrate_step, self.phi_range))

            return integrate_sum

        def calculate_L_nu_on_energy(self, energy):
            # КэВ
            ''' распределение L_nu от фазы на какой-то энергии излучения '''
            frequency = accretionColumnService.get_frequency_from_energy(energy)

            dS_simps = self.create_ds_for_integral()
            plank_func = accretionColumnService.plank_energy_on_frequency(frequency, self.T_eff)

            integrate_step = [0] * config.N_phi_accretion
            integrate_sum = [0] * config.t_max
            for rotation_index in range(config.t_max):
                for phi_index in range(config.N_phi_accretion):
                    integrate_step[phi_index] = 4 * np.pi * np.abs(scipy.integrate.simps(
                        plank_func * np.array(dS_simps) * np.array(
                            self.cos_psi_range[rotation_index][phi_index][:]), self.theta_range))
                integrate_sum[rotation_index] = np.abs(scipy.integrate.simps(integrate_step, self.phi_range))

            return integrate_sum

        def async_calculate_L_nu_on_energy(self, t_index, energy):
            frequency = accretionColumnService.get_frequency_from_energy(energy)

            dS_simps = self.create_ds_for_integral()
            plank_func = accretionColumnService.plank_energy_on_frequency(frequency, self.T_eff)

            integrate_step = [0] * config.N_phi_accretion

            for phi_index in range(config.N_phi_accretion):
                integrate_step[phi_index] = 4 * np.pi * np.abs(scipy.integrate.simps(
                    plank_func * np.array(dS_simps) * np.array(
                        self.cos_psi_range[t_index][phi_index][:]), self.theta_range))
            integrate_sum = np.abs(scipy.integrate.simps(integrate_step, self.phi_range))

            return integrate_sum

        def calculate_nu_L_nu_on_energy(self, energy):
            # КэВ
            frequency = accretionColumnService.get_frequency_from_energy(energy)
            return np.array(self.calculate_L_nu_on_energy(energy)) * frequency

        def calculate_black_body_approximation(self, energy, T_eff):
            # T_eff = np.mean(self.T_eff)
            frequency = accretionColumnService.get_frequency_from_energy(energy)

            dS_simps = self.create_ds_for_integral()
            plank_func = accretionColumnService.plank_energy_on_frequency(frequency, T_eff)

            integrate_step = [0] * config.N_phi_accretion
            integrate_sum = [0] * config.t_max
            for rotation_index in range(config.t_max):
                for phi_index in range(config.N_phi_accretion):
                    integrate_step[phi_index] = 4 * np.pi * np.abs(scipy.integrate.simps(
                        plank_func * np.array(dS_simps) * np.array(
                            self.cos_psi_range[rotation_index][phi_index][:]), self.theta_range))
                integrate_sum[rotation_index] = np.abs(scipy.integrate.simps(integrate_step, self.phi_range))

            return integrate_sum
