import numpy as np

import config


def get_angles_from_vector(vector):
    x = vector[0, 0]
    y = vector[0, 1]
    z = vector[0, 2]
    theta = np.arccos(z)  # np.arccos(z/r)
    if x > 0:
        if y >= 0:
            phi = np.arctan(y / x)
        else:
            phi = np.arctan(y / x) + 2 * np.pi
    elif x == 0:
        if y > 0:
            phi = np.pi / 2
        else:
            phi = 3 * np.pi / 2
    else:
        phi = np.arctan(y / x) + np.pi
    return phi, theta


def get_angles_2d(vector):
    x = vector[0, 0]
    y = vector[0, 1]
    z = vector[0, 2]
    r = (x ** 2 + y ** 2 + z ** 2) ** (1 / 2)
    theta = np.arccos(z / r)
    phi = np.arctan2(y, x)
    return phi, theta


def get_angles(vector):
    x = vector[0]
    y = vector[1]
    z = vector[2]
    r = (x ** 2 + y ** 2 + z ** 2) ** (1 / 2)
    theta = np.arccos(z / r)
    phi = np.arctan2(y, x)
    return phi, theta


def get_angles_from_vector_one_dimension(vector):
    x = vector[0]
    y = vector[1]
    z = vector[2]
    r = x ** 2 + y ** 2 + z ** 2
    theta = np.arccos(z / r)  # np.arccos(z/r)
    if x > 0:
        if y >= 0:
            phi = np.arctan(y / x)
        else:
            phi = np.arctan(y / x) + 2 * np.pi
    elif x == 0:
        if y > 0:
            phi = np.pi / 2
        else:
            phi = 3 * np.pi / 2
    else:
        phi = np.arctan(y / x) + np.pi
    return phi, theta


def get_angles_from_vector_one_dimension_with_r_ns(vector):
    x = vector[0] / config.R_ns
    y = vector[1] / config.R_ns
    z = vector[2] / config.R_ns
    theta = np.arccos(z)  # np.arccos(z/r)
    if x > 0:
        if y >= 0:
            phi = np.arctan(y / x)
        else:
            phi = np.arctan(y / x) + 2 * np.pi
    else:
        phi = np.arctan(y / x) + np.pi
    return phi, theta
