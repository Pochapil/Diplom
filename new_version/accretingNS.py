import numpy as np
import scipy.integrate
from scipy import interpolate

import newService
import config
import BS_approach as get_T_eff
import geometricTask.matrix as matrix

# from accretionColumn import AccretionColumn

column_surf_types = {'bot': 'bot', 'top': 'top'}
surface_surf_types = {'outer': 'outer', 'inner': 'inner'}


class AccretingPulsarConfiguration:

    def __init__(self, mu, beta_mu, mc2, a_portion, phi_0):
        self.mu = mu
        self.beta_mu = beta_mu
        self.mc2 = mc2
        self.a_portion = a_portion
        self.phi_0 = phi_0

        M_accretion_rate = mc2 * config.L_edd / config.c ** 2
        self.R_alfven = (mu ** 2 / (2 * M_accretion_rate * (2 * config.G * config.M_ns) ** (1 / 2))) ** (2 / 7)
        self.R_e = config.ksi_param * self.R_alfven

        self.top_column = AccretionColumn(self.R_e, mu, mc2, self.a_portion, self.phi_0,
                                          column_type=column_surf_types['top'])
        # Accret_col : equatorial, polar surfaces
        self.bot_column = AccretionColumn(self.R_e, mu, mc2, self.a_portion, self.phi_0,
                                          column_type=column_surf_types['bot'])

        self.top_magnet_lines = ...
        self.bot_magnet_lines = ...

        pass


class AccretionColumn:

    def __init__(self, R_e, mu, mc2, a_portion, phi_0, column_type):
        self.R_e = R_e
        self.column_type = column_type
        self.a_portion = a_portion
        self.phi_0 = phi_0

        # попробую что
        # R_e_outer_surface, R_e_inner_surface = self.R_e, self.R_e
        self.R_e_outer_surface, self.R_e_inner_surface = (1 + config.dRe_div_Re) * R_e, R_e  # допущение что толщина = 0

        M_accretion_rate = mc2 * config.L_edd / config.c ** 2
        self.T_eff, self.ksi_shock, self.L_x, self.beta = get_T_eff.get_Teff_distribution(
            newService.get_delta_distance_at_surface_NS(self.R_e),
            newService.get_A_normal_at_surface_NS(self.R_e, self.a_portion),
            mu, M_accretion_rate
        )

        self.outer_surface = Surface(self.R_e_outer_surface, self.a_portion, self.phi_0, self.ksi_shock,
                                     self.column_type, surface_type=surface_surf_types['outer'])
        self.inner_surface = Surface(self.R_e_inner_surface, self.a_portion, self.phi_0, self.ksi_shock,
                                     self.column_type, surface_type=surface_surf_types['inner'])

        self.correct_T_eff()

    def correct_T_eff(self):
        # так как нахожу распределение Teff по отрезку кси, то нужно привести к виду по сетке theta !!
        # нужно проинтерполировать внутри отрезка

        # создание сетки по которой было найдено распределение T_eff
        ksi_range = np.linspace(1, self.ksi_shock, config.N_theta_accretion)
        if self.ksi_shock > self.R_e / config.R_ns:
            ksi_range = np.linspace(1, self.R_e / config.R_ns, config.N_theta_accretion)

        # интерполяция
        x = ksi_range
        y = self.T_eff
        f = interpolate.interp1d(x, y, kind='cubic')

        # новый диапазон по х - начинаем с 1 индекса так как 0 элемент совпадает - x[0] = 1
        x_new = self.R_e / config.R_ns * np.sin(self.inner_surface.theta_range[1:-1]) ** 2
        y_new = f(x_new)

        # выравнивание
        for i in range(0, (len(y_new))):
            self.T_eff[i + 1] = y_new[i]


class Surface:
    def __init__(self, surf_R_e, a_portion, phi_0, ksi_shock, column_type, surface_type):
        self.surf_R_e = surf_R_e  # R_e на котором сидит - R_e либо R_e + \delta R_e
        self.surface_type = surface_type

        self.theta_accretion_begin = newService.get_theta_accretion_begin(self.surf_R_e)
        if (config.R_ns * ksi_shock / self.surf_R_e) >= 1:
            # есть набор параметров при которых модель не работает и ударная волна дальше магнитосферы, берем 90
            '''вопрос - мб нужна формула с arctan...'''
            self.theta_accretion_end = np.pi / 2
        else:
            # из усл силовой линии МП : r = R_e sin**2; end: ksi_shock = R_e sin**2
            self.theta_accretion_end = np.arcsin((config.R_ns * ksi_shock / self.surf_R_e) ** (1 / 2))

        phi_delta = 0
        # для нижней сместить углы
        if column_type == column_surf_types['bot']:
            self.theta_accretion_begin = np.pi - self.theta_accretion_begin
            self.theta_accretion_end = np.pi - self.theta_accretion_end
            phi_delta = np.pi

        self.theta_range = np.linspace(self.theta_accretion_begin, self.theta_accretion_end, config.N_theta_accretion)
        # phi_0 - центр колонки!!!!!!!!!!!!
        self.phi_range = np.linspace(-np.pi * a_portion, np.pi * a_portion, config.N_phi_accretion) + phi_0 + phi_delta
        if config.FLAG_PHI_0_OLD:
            # - в терминах старого phi
            self.phi_range = np.linspace(0, 2 * np.pi * a_portion) + phi_0 + phi_delta

        self.array_normal = self.create_array_normal(self.phi_range, self.theta_range, self.surface_type)

    def create_array_normal(self, phi_range, theta_range, surface_type):
        # array_normal = []  # матрица нормалей
        array_normal = np.zeros((config.N_phi_accretion, config.N_theta_accretion), dtype=object)
        coefficient = -1
        if surface_type == surface_surf_types['outer']:  # True - внешняя поверхность, False - внутренняя
            coefficient = 1
        for i in range(config.N_phi_accretion):
            for j in range(config.N_theta_accretion):
                # array_normal.append(coefficient * matrix.newE_n(phi_range[i], theta_range[j]))
                array_normal[i, j] = coefficient * matrix.newE_n(phi_range[i], theta_range[j])
        return array_normal

    def create_ds_for_integral(self):
        pass


class MagnetLine:
    def __init__(self):
        ...
