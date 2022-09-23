import numpy as np

cos = np.dot([1, 0, 1], [0, 1, 1])
print(cos)

from scipy.optimize import fsolve

R_e = 2

x_origin = 0
y_origin = 1
z_origin = 0

x_direction = 1
y_direction = 1
z_direction = 1


def equations(variables):
    theta, phi, t = variables
    eqn = np.empty(3)

    x = x_origin + x_direction * t
    y = y_origin + y_direction * t
    z = z_origin + z_direction * t

    eqn[0] = R_e * np.sin(theta) ** 3 * np.cos(phi) - x
    eqn[1] = R_e * np.sin(theta) ** 3 * np.sin(phi) - y
    eqn[2] = R_e * np.sin(theta) ** 2 * np.cos(theta) - z

    return eqn


# result_theta, result_phi, result_t = fsolve(equations, (0, 0, 0))
result = fsolve(equations, (0, 0, 0))
for i in result:
    print(i)
# print(result_theta, result_phi, result_t)

print(equations(result))

print(x_origin + x_direction * result[2])
print(y_origin + y_direction * result[2])
print(z_origin + z_direction * result[2])


# def equations(variables):
#     theta, phi, t = variables
#     eqn = np.empty(3)
#
#     x = x_origin + x_direction * t
#     y = y_origin + y_direction * t
#     z = z_origin + z_direction * t
#
#     eqn[0] = R_e * np.sin(theta) ** 3 * np.cos(phi) - x
#     eqn[1] = R_e * np.sin(theta) ** 3 * np.sin(phi) - y
#     eqn[2] = R_e * np.sin(theta) ** 2 * np.cos(theta) - z
#
#     return eqn
#
#     # result_theta, result_phi, result_t = fsolve(equations, (0, 0, 0))
#
#
# result = fsolve(equations, (origin_theta, origin_phi, 0))
#
# result_theta = result[0] - np.floor(result[0] / np.pi / 2) * 2 * np.pi
# result_phi = result[1] - np.floor(result[1] / np.pi / 2) * 2 * np.pi
#
# if result[2] > 0 and (
#         theta_accretion_begin < result_theta < theta_accretion_end or theta_accretion_end_1 < result_theta < theta_accretion_begin_1) and result_phi < lim_phi_accretion:
#     return True


'''new main'''
# для не 0 угла наблюдателя по фи
        # eta = np.sin(direction_theta) / np.sin(origin_theta)
        # cos_alpha = np.sin(origin_theta) * np.cos(origin_phi) * np.sin(direction_theta) * np.cos(direction_phi) \
        #             + np.sin(origin_theta) * np.sin(origin_phi) * np.sin(direction_theta) * np.sin(direction_phi) \
        #             + np.cos(origin_theta) * np.cos(direction_theta)
        #
        # c_x_5 = 1
        # c_x_4 = 6 * cos_alpha
        # c_x_3 = 3 + 12 * cos_alpha ** 2 - eta ** 4
        # c_x_2 = 12 * cos_alpha + 8 * cos_alpha ** 3 - 4 * np.cos(origin_phi) * np.cos(direction_phi) * eta ** 3 \
        #         - 4 * np.sin(origin_phi) * np.sin(direction_phi) * eta ** 3
        # c_x_1 = 3 + 12 * cos_alpha ** 2 \
        #         - 2 * eta ** 2 \
        #         - 4 * np.cos(origin_phi) ** 2 * np.cos(direction_phi) ** 2 * eta ** 2 \
        #         - 4 * np.sin(origin_phi) ** 2 * np.sin(direction_phi) ** 2 * eta ** 2 \
        #         - 8 * np.cos(origin_phi) * np.cos(direction_phi) * np.sin(origin_phi) * np.sin(direction_phi) * eta ** 2
        # c_x_0 = 6 * cos_alpha - 4 * np.cos(origin_phi) * np.cos(direction_phi) * eta \
        #         - 4 * np.sin(origin_phi) * np.sin(direction_phi) * eta


# r = R_e * np.sin(origin_theta) ** 2
#         # декартовая СК из сферических
#         origin_x = np.sin(origin_theta) * np.cos(origin_phi) * r
#         origin_y = np.sin(origin_theta) * np.sin(origin_phi) * r
#         origin_z = np.cos(origin_theta) * r
# так как находили корни в другой СК то нужно перевести все снова
        # observer_coord_sys_direction_x = 1 * np.sin(direction_theta)
        # observer_coord_sys_direction_y = 0
        # observer_coord_sys_direction_z = 1 * np.cos(direction_theta)
        #
        # observer_coord_sys_origin_x = np.sin(origin_theta) * np.cos(phi_delta) * r
        # observer_coord_sys_origin_y = np.sin(origin_theta) * np.sin(phi_delta) * r
        # observer_coord_sys_origin_z = np.cos(origin_theta) * r
        #
        # direction_x = observer_coord_sys_direction_x
        # direction_y = observer_coord_sys_direction_y
        # direction_z = observer_coord_sys_direction_z
        #
        # origin_x = observer_coord_sys_origin_x
        # origin_y = observer_coord_sys_origin_y
        # origin_z = observer_coord_sys_origin_z