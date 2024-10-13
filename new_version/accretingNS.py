import numpy as np
import scipy.integrate
from scipy import interpolate

import accretionColumnService
import config
import BS_distribution_T_eff as get_T_eff
import geometricTask.matrix as matrix


# from accretionColumn import AccretionColumn


class AccretionColumn:

    def __init__(self, R_e, column_type):
        self.R_e_outer_surface, self.R_e_inner_surface = R_e, R_e  # допущение что толщина = 0
        # попробую что
        # R_e_outer_surface, R_e_inner_surface = self.R_e, self.R_e

        self.theta_accretion_begin_outer_surface = accretionColumnService.get_theta_accretion_begin(
            self.R_e_outer_surface)
        self.theta_accretion_begin_inner_surface = accretionColumnService.get_theta_accretion_begin(
            self.R_e_inner_surface)

        self.T_eff, self.ksi_shock, self.L_x, self.beta = get_T_eff.get_Teff_distribution()

        if column_type:
            self.theta_accretion_begin_outer_surface = np.pi - self.theta_accretion_begin_outer_surface
            self.theta_accretion_begin_inner_surface = np.pi - self.theta_accretion_begin_inner_surface

        self.outer_surface = Surface(self.theta_accretion_begin_outer_surface, self.R_e_outer_surface, True,
                                     self.column_type)
        self.inner_surface = Surface(self.theta_accretion_begin_inner_surface, self.R_e_inner_surface, False,
                                     self.column_type)


class Surface:
    def __init__(self, theta_accretion_begin, surf_R_e, surface_type, column_type):
        self.array_normal = self.create_array_normal()
        self.R_magnet_line = surf_R_e  # R_e на котором сидит - R_e либо R_e + \delta R_e

        phi_delta = 0

        if (config.R_ns * self.ksi_shock / self.R_magnet_line) >= 1:
            # есть набор параметров при которых модель не работает и ударная волна дальше магнитосферы, берем 90
            '''вопрос - мб нужна формула с arctan...'''
            theta_accretion_end = np.pi / 2
        else:
            # из усл силовой линии МП : r = R_e sin**2; end: ksi_shock = R_e sin**2
            theta_accretion_end = np.arcsin((config.R_ns * self.ksi_shock / self.R_magnet_line) ** (1 / 2))

        if not column_type:  # column_type: True - top, False - bot
            # для нижней сместить углы
            # theta_begin было смещено в main_loop 116 строка!!
            theta_accretion_end = np.pi - theta_accretion_end
            phi_delta = np.pi

    def create_array_normal(self, phi_range, theta_range, surface_type=True):
        # array_normal = []  # матрица нормалей
        array_normal = np.zeros((config.N_phi_accretion, config.N_theta_accretion), dtype=object)
        coefficient = -1
        if surface_type:  # True - внешняя поверхность, False - внутренняя
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


class AccretingPulsarConfiguration:

    def __init__(self, mu, beta, M_accretion_rate, a, phi_0, R_e_outer_surface, theta_accretion_begin_outer_surface,
                 R_e_inner_surface, theta_accretion_begin_inner_surface):
        self.R_alfven = (mu ** 2 / (2 * M_accretion_rate * (2 * config.G * config.M_ns) ** (1 / 2))) ** (2 / 7)
        self.R_e = config.ksi_param * self.R_alfven

        self.top_column = AccretionColumn(column_type=1)
        # Accret_col : equatorial, polar surfaces
        self.bot_column = AccretionColumn(column_type=-1)

        self.top_magnet_lines = ...
        self.bot_magnet_lines = ...

        pass
