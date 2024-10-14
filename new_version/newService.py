import numpy as np

import vectors
import config

approx_type = True
approx_method = "cone"
if approx_type:
    approx_method = "dipole"


def get_delta_distance_at_surface_NS(R_e):
    # R=R_e * sin_theta ** 2
    theta = get_theta_accretion_begin(R_e)
    print(theta)
    return get_delta_distance(theta, R_e)


def get_A_normal_at_surface_NS(R_e, a_portion):
    theta = get_theta_accretion_begin(R_e)
    print(theta)
    # a - в азимутальном направлении поток занимает фиксированную долю a полного круга 2πR sinθ
    return get_A_normal(theta, R_e, a_portion)


# формула 2 в статье
# the distance between the field lines δ
# ширина аккреции на поверхности м
def get_delta_distance(theta, R_e):
    # R=R_e * sin_theta ** 2
    return R_e * np.sin(theta) ** 3 / (1 + 3 * np.cos(theta) ** 2) ** (1 / 2) * config.dRe_div_Re


# формула 3 в статье
# площадь поперечного сечения аккреционного потока
# a cross-section
def get_A_normal(theta, R_e, a_portion):
    # a - в азимутальном направлении поток занимает фиксированную долю a полного круга 2πR sinθ
    return 2 * get_delta_distance(theta, R_e) * 2 * np.pi * a_portion * R_e * np.sin(theta) ** 3


# из усл силовой линии МП : r = R_e sin**2
def get_theta_accretion_begin(R_e):
    return np.arcsin((config.R_ns / R_e) ** (1 / 2))


def find_intersect_solution(a, b, c):
    if b ** 2 - 4 * a * c >= 0:
        t_1 = (-b + (b ** 2 - 4 * a * c) ** (1 / 2)) / (2 * a)
        t_2 = (-b - (b ** 2 - 4 * a * c) ** (1 / 2)) / (2 * a)
        return t_1, t_2
    else:
        return -1, -1


def get_free_fall_velocity(theta, R_e):
    r = R_e * np.sin(theta) ** 2
    return (2 * config.G * config.M_ns / r) ** (1 / 2)


def get_tau_for_opacity(theta, R_e, M_accretion_rate, a_portion):
    tau = config.k * M_accretion_rate * get_delta_distance(theta, R_e) / (
            get_A_normal(theta, R_e, a_portion) * get_free_fall_velocity(theta, R_e))
    return tau


def intersection_with_sphere(origin_x, origin_y, origin_z, direction_x, direction_y, direction_z):
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


def intersection_with_cone(origin_x, origin_y, origin_z, direction_x, direction_y, direction_z, ksi_shock,
                           theta_accretion_begin, theta_accretion_end, top_column_phi_range, bot_column_phi_range,
                           surface_type):
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


def get_pulsed_fraction(arr):
    min_value = min(arr)
    max_value = max(arr)
    if max_value == 0:
        return 0
    return (max_value - min_value) / (max_value + min_value)


# нужно ли на пи?
def plank_energy_on_wavelength(wavelength, T):
    return 2 * config.h_plank_ergs * config.c ** 2 / wavelength ** 5 \
           * 1 / (np.e ** (config.h_plank_ergs * config.c / (wavelength * config.k_bolc * T)) - 1)


def plank_energy_on_frequency(frequency, T):
    # erg / s / sr / cm**2 / hz
    return 2 * config.h_plank_ergs * frequency ** 3 / config.c ** 2 \
           * 1 / (np.e ** (config.h_plank_ergs * frequency / (config.k_bolc * T)) - 1)


def get_frequency_from_energy(energy):
    coefficient = 1000  # КэВ а не эВ
    frequency = coefficient * energy / config.h_plank_evs  # E = h f
    return frequency  # гц!
