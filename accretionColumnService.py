import numpy as np

import config
import vectors

approx_type = True
approx_method = "cone"
if approx_type:
    approx_method = "dipole"


# формула 2 в статье
# the distance between the field lines δ
# ширина аккреции на поверхности м
def get_delta_distance(theta, R_e):
    # R=R_e * sin_theta ** 2
    return R_e * np.sin(theta) ** 3 / (1 + 3 * np.cos(theta) ** 2) ** (1 / 2) * config.dRe_div_Re


# формула 3 в статье
# площадь поперечного сечения аккреционного потока
# a cross-section
def get_A_normal(theta, R_e):
    # a - в азимутальном направлении поток занимает фиксированную долю a полного круга 2πR sinθ
    return 2 * get_delta_distance(theta, R_e) * 2 * np.pi * config.a_portion * R_e * np.sin(theta) ** 3


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


def get_tau_for_opacity(theta, R_e):
    tau = config.k * config.M_accretion_rate * get_delta_distance(theta, R_e) / (
            get_A_normal(theta, R_e) * get_free_fall_velocity(theta, R_e))
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


def intersection_with_dipole_lines(origin_phi, origin_theta, direction_vector, origin_x, origin_y, origin_z,
                                   direction_x, direction_y, direction_z, ksi_shock, theta_accretion_end,
                                   top_column_phi_range, bot_column_phi_range, r):
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

    direction_phi, direction_theta = vectors.get_angles_from_vector(direction_vector)

    phi_delta = origin_phi - direction_phi
    eta = np.sin(direction_theta) / np.sin(origin_theta)
    cos_alpha = np.sin(origin_theta) * np.cos(phi_delta) * np.sin(direction_theta) + np.cos(origin_theta) * np.cos(
        direction_theta)

    '''коэффициенты получены аналитически приравнивая уравнение луча и уравнение дипольной линии (ур. 4.46 в дипломе)'''
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
            # if not intersect_r_correct:
            #     print('intersect_r not correct')
            # для верхней колонки:
            if theta_accretion_end > np.arcsin((2 / 3) ** (1 / 2)):
                R_alfven = (config.mu ** 2 / (
                        2 * config.M_accretion_rate * (2 * config.G * config.M_ns) ** (1 / 2))) ** (2 / 7)
                R_e = config.ksi_param * R_alfven
                intersect_z_correct = 0 <= intersect_point[2] <= R_e / config.R_ns * 0.3849001794597505097
            else:
                intersect_z_correct = 0 <= intersect_point[2] <= ksi_shock * np.cos(theta_accretion_end)
            intersect_phi_correct = (top_column_phi_range[0] <= intersect_phi <= top_column_phi_range[-1]) or (
                    0 <= intersect_phi <= top_column_phi_range[-1] - 2 * np.pi)
            if intersect_z_correct and intersect_phi_correct and intersect_r_correct:
                return True

            # для нижней колонки:
            if theta_accretion_end > np.arcsin((2 / 3) ** (1 / 2)):
                R_alfven = (config.mu ** 2 / (
                        2 * config.M_accretion_rate * (2 * config.G * config.M_ns) ** (1 / 2))) ** (2 / 7)
                R_e = config.ksi_param * R_alfven
                intersect_z_correct = - R_e / config.R_ns * 0.3849001794597505097 <= intersect_point[2] <= 0
            else:
                intersect_z_correct = - (ksi_shock * np.cos(theta_accretion_end)) <= intersect_point[2] <= 0
            # условие - так как углы из метода get_angles_from_vector_one_dimension от 0 до 2 pi поэтому
            # нужно учесть углы которые превышают 2 pi в 1 скобке условия
            intersect_phi_correct = (0 <= intersect_phi <= bot_column_phi_range[-1] - 2 * np.pi) or (
                    bot_column_phi_range[0] <= intersect_phi <= bot_column_phi_range[-1])
            if intersect_z_correct and intersect_phi_correct and intersect_r_correct:
                return True

    return False


def check_intersection_with_sphere_and_columns(origin_phi, origin_theta, direction_vector, ksi_shock,
                                               theta_accretion_begin, theta_accretion_end, top_column_phi_range,
                                               bot_column_phi_range, surface_type, R_e):
    # сначала поиск со сферой после - с конусом (а нужно ли со сферой искать ?)
    # все нормирую на радиус НЗ
    # r = origin + t * direction - уравнение луча

    if top_column_phi_range[0] - 2 * np.pi > 0:
        top_column_phi_range = top_column_phi_range - 2 * np.pi

    if bot_column_phi_range[0] - 2 * np.pi > 0:
        bot_column_phi_range = bot_column_phi_range - 2 * np.pi

    r = R_e / config.R_ns * np.sin(origin_theta) ** 2
    # декартовая СК из сферических
    origin_x = np.sin(origin_theta) * np.cos(origin_phi) * r
    origin_y = np.sin(origin_theta) * np.sin(origin_phi) * r
    origin_z = np.cos(origin_theta) * r
    # origin_point = [x, y, z]

    direction_x = direction_vector[0, 0]
    direction_y = direction_vector[0, 1]
    direction_z = direction_vector[0, 2]

    # direction_phi, direction_theta = vectors.get_angles_from_vector(direction_vector)

    if approx_type:
        return (intersection_with_sphere(origin_x, origin_y, origin_z, direction_x, direction_y, direction_z) or
                intersection_with_dipole_lines(origin_phi, origin_theta, direction_vector, origin_x, origin_y, origin_z,
                                               direction_x, direction_y, direction_z, ksi_shock, theta_accretion_end,
                                               top_column_phi_range, bot_column_phi_range, r))
    else:
        return (intersection_with_sphere(origin_x, origin_y, origin_z, direction_x, direction_y, direction_z) or
                intersection_with_cone(origin_x, origin_y, origin_z, direction_x, direction_y, direction_z, ksi_shock,
                                       theta_accretion_begin, theta_accretion_end, top_column_phi_range,
                                       bot_column_phi_range,
                                       surface_type))


def intersection_with_dipole_lines_with_opacity(origin_phi, origin_theta, direction_vector, origin_x, origin_y,
                                                origin_z, direction_x, direction_y, direction_z, ksi_shock,
                                                theta_accretion_end, top_column_phi_range, bot_column_phi_range, R_e,
                                                r):
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

    direction_phi, direction_theta = vectors.get_angles_from_vector(direction_vector)

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

            intersect_r_correct = intersect_r < R_e / config.R_ns

            # для верхней колонки:
            # 0.3849 - max z: (sin**2 * cos == 0) => x = arcsin((2/3)**(1/2)) => sin**2 x * cos x = 0.3849
            intersect_z_correct = 0 < intersect_point[2] <= R_e / config.R_ns * 0.3849001794597505097
            intersect_phi_correct = (top_column_phi_range[0] <= intersect_phi <= top_column_phi_range[-1]) or (
                    0 <= intersect_phi <= top_column_phi_range[-1] - 2 * np.pi)
            if intersect_z_correct and intersect_phi_correct and intersect_r_correct:
                return True

            # для нижней колонки:
            intersect_z_correct = -(R_e / config.R_ns * 0.3849001794597505097) <= intersect_point[2] < 0
            # условие - так как углы из метода get_angles_from_vector_one_dimension от 0 до 2 pi поэтому
            # нужно учесть углы которые превышают 2 pi в 1 скобке условия
            intersect_phi_correct = (0 <= intersect_phi <= bot_column_phi_range[-1] - 2 * np.pi) or (
                    bot_column_phi_range[0] <= intersect_phi <= bot_column_phi_range[-1])
            if intersect_z_correct and intersect_phi_correct and intersect_r_correct:
                return True

    return False


def get_solutions_for_dipole_magnet_lines(origin_phi, origin_theta, direction_vector):
    # для 0 угла наблюдателя по фи. поэтому находим phi_delta
    direction_phi, direction_theta = vectors.get_angles_2d(direction_vector)
    # get_angles_from_vector(direction_vector)
    if direction_phi < 0:
        direction_phi += 2 * np.pi
    # вспомогательные переменные, были введены для упрощения аналитического вывода
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

    return solutions


def get_vals(origin_phi, origin_theta, direction_vector, origin_x, origin_y,
             origin_z, direction_x, direction_y, direction_z, ksi_shock,
             top_column_phi_range, bot_column_phi_range, top_column_theta_range, bot_column_theta_range,
             R_e, r, betta_mu):
    '''
    есть аналитическое уравнение для полинома дипольной линии 5 степени
    находим уравнение в сферических координатах.

    достаем корни
    ищем положительные
    находим пересечение с колонками
        если пересекается то истина
        если нет то ложь

    корни в магнитной СК - betta_mu
    '''

    solutions = get_solutions_for_dipole_magnet_lines(origin_phi, origin_theta, direction_vector)

    '''думаю что нужно только phi theta для нахождения правильных пересечений. r,z нужно было раньше когда проверял на 
    пересечение с колонкой'''
    # мб добавить условие на минимум t.real -- and solution.real > 1e-2
    for solution in solutions:
        if solution.real > 0 and solution.imag == 0:
            direction_t = solution.real * r
            intersect_point = np.array([origin_x, origin_y, origin_z]) + direction_t * np.array(
                [direction_x, direction_y, direction_z])

            intersect_phi, intersect_theta = vectors.get_angles(intersect_point)
            if intersect_phi < 0:
                intersect_phi += 2 * np.pi

            intersect_r = (intersect_point[0] ** 2 + intersect_point[1] ** 2 + intersect_point[2] ** 2) ** (1 / 2)

            # if not (abs(R_e / config.R_ns * np.sin(intersect_theta) ** 2 - intersect_r) < 0.001):
            #     print(abs(R_e / config.R_ns * np.sin(intersect_theta) ** 2 - intersect_r))

            # curr = abs(R_e / config.R_ns * np.sin(intersect_theta) ** 2 - intersect_r)
            # новые линии магнитосферы
            # theta_end = np.pi / 2 - betta_mu * np.cos(intersect_phi) - ОШИБКА !!
            theta_end = np.pi / 2 - np.arctan(np.tan(betta_mu) * np.cos(intersect_phi))

            # для верхней колонки:
            top_column_intersect_phi_correct = (top_column_phi_range[0] <= intersect_phi <= top_column_phi_range[
                -1]) or (0 <= intersect_phi <= top_column_phi_range[-1] - 2 * np.pi)
            top_column_intersect_theta_correct = intersect_theta < theta_end

            # для нижней колонки:
            bot_column_intersect_phi_correct = (bot_column_phi_range[0] <= intersect_phi <= bot_column_phi_range[
                -1]) or (0 <= intersect_phi <= bot_column_phi_range[-1] - 2 * np.pi)
            bot_column_intersect_theta_correct = intersect_theta > theta_end

            # проверяем на пересечение с колонками
            # intersect_r_correct = intersect_r > ksi_shock
            # if not intersect_r_correct and (top_column_intersect_phi_correct or bot_column_intersect_phi_correct):
            #     return 0

            if (intersect_theta < top_column_theta_range[-1] and top_column_intersect_phi_correct) or (
                    intersect_theta > bot_column_theta_range[-1] and bot_column_intersect_phi_correct):
                return 0

            intersection_condition = (top_column_intersect_phi_correct and top_column_intersect_theta_correct) or (
                    bot_column_intersect_phi_correct and bot_column_intersect_theta_correct)

            # if not (top_column_intersect_phi_correct and top_column_intersect_theta_correct) and (
            #         bot_column_intersect_phi_correct and bot_column_intersect_theta_correct):
            #     print('пересекли bot в точке')
            #     print(intersect_phi / config.grad_to_rad, intersect_theta / config.grad_to_rad)
            #     print(f'theta_end={theta_end / config.grad_to_rad}')

            if config.tau_flag:
                if intersection_condition:
                    tau = get_tau_for_opacity(intersect_theta, R_e)
                    if tau > config.tau_cutoff:
                        return np.exp(-1 * tau)
                    else:
                        return 1

            elif config.opacity_above_shock > 0:
                if intersection_condition:
                    return 1 - config.opacity_above_shock

            else:
                return 1
    return 1


def get_data_for_magnet_lines(theta_range_column, phi_range_column, fi_0):
    '''смог реализовать только с помощью маски'''
    theta_array_end = np.pi / 2 + config.betta_mu
    # ограничиваю колонкой
    theta_array_end = min((np.pi - theta_range_column[-1]), theta_array_end)
    theta_array_begin = theta_range_column[-1]

    step_theta_accretion = (theta_array_end - theta_array_begin) / (config.N_theta_accretion - 1)
    theta_range = np.array(
        [theta_array_begin + step_theta_accretion * j for j in range(config.N_theta_accretion)])

    step_phi_accretion = config.lim_phi_accretion / (config.N_phi_accretion - 1)
    phi_range = np.array(
        [fi_0 * config.grad_to_rad + step_phi_accretion * i for i in range(config.N_phi_accretion)])

    r, p = np.meshgrid(np.sin(theta_range) ** 2, phi_range)
    r1 = r * np.sin(theta_range)
    x = r1 * np.cos(p)
    y = r1 * np.sin(p)
    z = r * np.cos(theta_range)

    mask = np.zeros_like(x).astype(bool)
    for i in range(len(phi_range)):
        for j in range(len(theta_range)):
            theta_end = np.pi / 2 - config.betta_mu * np.cos(phi_range[i])
            if theta_range[j] > theta_end:
                # pass
                mask[i][j] = True

    return x, y, z, mask


def get_attenuation_coefficient_intersection_with_dipole_tau(origin_phi, origin_theta, direction_vector, origin_x,
                                                             origin_y, origin_z, direction_x, direction_y, direction_z,
                                                             ksi_shock, theta_accretion_end, top_column_phi_range,
                                                             bot_column_phi_range, R_e, r):
    direction_phi, direction_theta = vectors.get_angles_from_vector(direction_vector)

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

            intersect_r_correct = intersect_r < R_e / config.R_ns

            # для верхней колонки:
            # 0.3849 - max z: (sin**2 * cos == 0) => x = arcsin((2/3)**(1/2)) => sin**2 x * cos x = 0.3849
            intersect_z_correct = 0 < intersect_point[2] <= R_e / config.R_ns * 0.3849001794597505097
            intersect_phi_correct = (top_column_phi_range[0] <= intersect_phi <= top_column_phi_range[-1]) or (
                    0 <= intersect_phi <= top_column_phi_range[-1] - 2 * np.pi)
            if intersect_z_correct and intersect_phi_correct and intersect_r_correct:
                if intersect_r < ksi_shock:
                    return 0
                else:
                    tau = get_tau_for_opacity(intersect_theta, R_e)
                    if tau > config.tau_cutoff:
                        return np.exp(-1 * tau)
                    else:
                        return 1

            # для нижней колонки:
            intersect_z_correct = -(R_e / config.R_ns * 0.3849001794597505097) <= intersect_point[2] < 0
            # условие - так как углы из метода get_angles_from_vector_one_dimension от 0 до 2 pi поэтому
            # нужно учесть углы которые превышают 2 pi в 1 скобке условия
            intersect_phi_correct = (0 <= intersect_phi <= bot_column_phi_range[-1] - 2 * np.pi) or (
                    bot_column_phi_range[0] <= intersect_phi <= bot_column_phi_range[-1])
            if intersect_z_correct and intersect_phi_correct and intersect_r_correct:
                if intersect_r < ksi_shock:
                    return 0
                else:
                    tau = get_tau_for_opacity(intersect_theta, R_e)
                    if tau > config.tau_cutoff:
                        return np.exp(-1 * tau)
                    else:
                        return 1

    return 1


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
