import geometricTask.matrix as matrix
import numpy as np
import matplotlib.pyplot as plt

import scipy.integrate

import matplotlib.cm as cm
from matplotlib.colors import Normalize

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


# от поверхности NS - угол при котором радиус = радиусу НЗ
theta_accretion_begin_outer_surface = get_theta_accretion_begin(R_e_outer_surface)
theta_accretion_begin_inner_surface = get_theta_accretion_begin(R_e_inner_surface)


class AccretionColumn:

    def __init__(self, R_e_outer_surface, theta_accretion_begin_outer_surface, R_e_inner_surface,
                 theta_accretion_begin_inner_surface):
        self.outer_surface = self.Surface(theta_accretion_begin_outer_surface, R_e_outer_surface, True)
        self.inner_surface = self.Surface(theta_accretion_begin_inner_surface, R_e_inner_surface, False)

    class Surface:
        def __init__(self, theta_accretion_begin, R_e, flag):
            self.R_e = R_e
            self.flag = flag

            self.T_eff, self.ksi_shock, self.L_x = get_T_eff.get_Teff_distribution(config.N_theta_accretion, R_e,
                                                                                   get_delta_distance(
                                                                                       theta_accretion_begin),
                                                                                   get_A_normal(theta_accretion_begin))
            theta_accretion_end = np.arcsin((config.R_ns * self.ksi_shock / R_e) ** (1 / 2))

            step_phi_accretion = config.lim_phi_accretion / config.N_phi_accretion
            step_theta_accretion = (theta_accretion_end - theta_accretion_begin) / config.N_theta_accretion

            self.theta_range = np.array(
                [theta_accretion_begin + step_theta_accretion * j for j in range(config.N_theta_accretion)])
            self.phi_range = np.array(
                [config.phi_accretion_begin + step_phi_accretion * i for i in range(config.N_phi_accretion)])

            self.cos_psi_range = np.empty([config.N_phi_accretion, config.N_theta_accretion])
            '''тут создать матрицу косинусов 1 раз и использовать потом'''

            self.array_normal = self.create_array_normal(self.phi_range, self.theta_range, self.flag)

        def create_array_normal(self, phi_range, theta_range, flag=True):
            array_normal = []  # матрица нормалей
            coefficient = -1
            if flag:  # True - внешняя поверхность, False - внутренняя
                coefficient = 1
            # matrix_normal = np.empty([N_phi_accretion, N_theta_accretion])
            for i in range(config.N_phi_accretion):
                for j in range(config.N_theta_accretion):
                    # matrix_normal[i, j] = matrix.newE_n(phi_range[i], theta_range[j])
                    array_normal.append(coefficient * matrix.newE_n(phi_range[i], theta_range[j]))
            return array_normal

